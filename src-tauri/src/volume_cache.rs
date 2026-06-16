//! Bounded LRU cache for decoded CT volumes.
//!
//! Keyed by `(dicom_dir, series_uid)`; holds pipeline metadata, an `Arc<Vec<i16>>`
//! of raw voxels (so the framed IPC response can be rebuilt without redecoding),
//! and the fully-assembled `LoadedVolume` (f32 `Array3` shared via `Arc`) that
//! downstream CPR / FAI / MMD pipelines operate on.
//!
//! The cache is kept small (default 3 entries) so worst-case resident memory
//! scales like `3 * sizeof(volume) ≈ 450 MB i16 + 900 MB f32 ≈ 1.35 GB`.
//! The 3 GB cornerstone cache on the frontend is the larger line item.

use std::collections::{HashMap, VecDeque};
use std::sync::Arc;

use pcat_pipeline::dicom_load::VolumeMetadata as PipelineVolumeMetadata;
use pcat_pipeline::types::LoadedVolume;

/// Default maximum number of decoded volumes retained in the LRU cache.
///
/// Sized for a full patient load: CaScore + CCTA + up to four MonoPlus keV
/// series resident at once so the header volume switcher can cross-reference
/// modalities without redecoding.
pub const VOLUME_CACHE_MAX: usize = 8;

/// A single cached volume. `voxels_i16` is shared via `Arc` so both the
/// IPC framed-response path and the cache entry reference the same buffer;
/// cloning the `Arc` is a refcount bump, not a ~150 MB data copy.
pub struct CachedVolume {
    pub metadata: PipelineVolumeMetadata,
    pub voxels_i16: Arc<Vec<i16>>,
    pub volume: LoadedVolume,
}

/// Bounded LRU cache keyed by `(dicom_dir, series_uid)`.
pub struct VolumeCache {
    inner: HashMap<(String, String), Arc<CachedVolume>>,
    /// LRU order. Oldest key at the front, most-recently-used at the back.
    order: VecDeque<(String, String)>,
    max_entries: usize,
}

impl VolumeCache {
    pub fn new(max_entries: usize) -> Self {
        Self {
            inner: HashMap::new(),
            order: VecDeque::new(),
            max_entries,
        }
    }

    /// Lookup by key. On hit, moves the entry to the back of the LRU queue.
    pub fn get(&mut self, key: &(String, String)) -> Option<Arc<CachedVolume>> {
        let hit = self.inner.get(key)?.clone();
        // Touch: move this key to the back of the order deque.
        if let Some(pos) = self.order.iter().position(|k| k == key) {
            if let Some(k) = self.order.remove(pos) {
                self.order.push_back(k);
            }
        }
        Some(hit)
    }

    /// Insert or replace an entry. Evicts the oldest entry if capacity is
    /// exceeded. Re-inserting an existing key refreshes its LRU position.
    pub fn insert(&mut self, key: (String, String), value: CachedVolume) {
        if self.inner.contains_key(&key) {
            // Replace existing entry: drop from order first, then re-push.
            if let Some(pos) = self.order.iter().position(|k| k == &key) {
                self.order.remove(pos);
            }
        }
        self.inner.insert(key.clone(), Arc::new(value));
        self.order.push_back(key);

        // Evict oldest while over capacity.
        while self.order.len() > self.max_entries {
            if let Some(oldest) = self.order.pop_front() {
                self.inner.remove(&oldest);
            }
        }
    }

    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.inner.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array3;

    /// Build a minimal `CachedVolume` sufficient for cache-behavior tests.
    /// Does NOT need to represent a real DICOM volume — we only exercise key
    /// lookup / LRU order here.
    fn make_dummy_cached_volume() -> CachedVolume {
        let meta = PipelineVolumeMetadata {
            series_uid: String::new(),
            series_description: String::new(),
            image_comments: None,
            rows: 1,
            cols: 1,
            num_slices: 1,
            pixel_spacing: [1.0, 1.0],
            slice_spacing: 1.0,
            orientation: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            window_center: 0.0,
            window_width: 1.0,
            patient_name: String::new(),
            study_description: String::new(),
            slice_positions_z: vec![0.0],
            image_position_patient: [0.0, 0.0, 0.0],
        };
        let arr = Array3::<f32>::zeros((1, 1, 1));
        let volume = LoadedVolume {
            data: Arc::new(arr),
            spacing: [1.0, 1.0, 1.0],
            origin: [0.0, 0.0, 0.0],
            direction: [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            window_center: 0.0,
            window_width: 1.0,
            patient_name: String::new(),
            study_description: String::new(),
        };
        CachedVolume {
            metadata: meta,
            voxels_i16: Arc::new(vec![0i16]),
            volume,
        }
    }

    #[test]
    fn cache_get_miss_returns_none() {
        let mut c = VolumeCache::new(3);
        assert!(c.get(&("a".into(), "x".into())).is_none());
    }

    #[test]
    fn cache_insert_then_get_hits() {
        let mut c = VolumeCache::new(3);
        c.insert(("a".into(), "x".into()), make_dummy_cached_volume());
        assert!(c.get(&("a".into(), "x".into())).is_some());
    }

    #[test]
    fn cache_lru_evicts_oldest_past_max() {
        let mut c = VolumeCache::new(2);
        c.insert(("a".into(), "x".into()), make_dummy_cached_volume());
        c.insert(("b".into(), "y".into()), make_dummy_cached_volume());
        c.insert(("c".into(), "z".into()), make_dummy_cached_volume());
        assert!(c.get(&("a".into(), "x".into())).is_none(), "a should be evicted");
        assert!(c.get(&("b".into(), "y".into())).is_some());
        assert!(c.get(&("c".into(), "z".into())).is_some());
    }

    #[test]
    fn cache_get_touches_lru_order() {
        // After inserting a, b and then getting a, inserting c should evict b, not a.
        let mut c = VolumeCache::new(2);
        c.insert(("a".into(), "x".into()), make_dummy_cached_volume());
        c.insert(("b".into(), "y".into()), make_dummy_cached_volume());
        let _ = c.get(&("a".into(), "x".into())); // touches a
        c.insert(("c".into(), "z".into()), make_dummy_cached_volume());
        assert!(c.get(&("a".into(), "x".into())).is_some(), "a was touched, should survive");
        assert!(c.get(&("b".into(), "y".into())).is_none(), "b should be evicted");
    }

    #[test]
    fn cache_reinsert_same_key_does_not_evict() {
        let mut c = VolumeCache::new(2);
        c.insert(("a".into(), "x".into()), make_dummy_cached_volume());
        c.insert(("b".into(), "y".into()), make_dummy_cached_volume());
        // Re-insert "a" — should refresh, not push a third entry.
        c.insert(("a".into(), "x".into()), make_dummy_cached_volume());
        assert_eq!(c.len(), 2);
        assert!(c.get(&("a".into(), "x".into())).is_some());
        assert!(c.get(&("b".into(), "y".into())).is_some());
    }
}
