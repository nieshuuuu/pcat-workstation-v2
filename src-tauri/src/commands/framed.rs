//! Framed binary response for Tauri commands.
//!
//! Layout (single buffer):
//!   [u32 LE: metadata_json_length] [metadata_json_bytes] [payload_bytes]
//!
//! Frontend parses the first 4 bytes, reads JSON, treats the rest as the
//! payload ArrayBuffer.

use serde::Serialize;

/// Encode (metadata, payload_bytes) into a single framed Vec<u8>.
pub fn encode_frame<M: Serialize>(metadata: &M, payload: &[u8]) -> Result<Vec<u8>, String> {
    let json = serde_json::to_vec(metadata).map_err(|e| format!("serialize metadata: {e}"))?;
    let meta_len: u32 = json
        .len()
        .try_into()
        .map_err(|_| "metadata json exceeds u32 range".to_string())?;

    let mut out = Vec::with_capacity(4 + json.len() + payload.len());
    out.extend_from_slice(&meta_len.to_le_bytes());
    out.extend_from_slice(&json);
    out.extend_from_slice(payload);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Serialize)]
    struct M { a: u32, b: String }

    #[test]
    fn round_trips_payload() {
        let m = M { a: 42, b: "hi".into() };
        let payload = [0x01u8, 0x02, 0x03, 0x04];
        let buf = encode_frame(&m, &payload).unwrap();

        let meta_len = u32::from_le_bytes(buf[..4].try_into().unwrap()) as usize;
        let json = &buf[4..4 + meta_len];
        let body = &buf[4 + meta_len..];
        let parsed: serde_json::Value = serde_json::from_slice(json).unwrap();
        assert_eq!(parsed["a"], 42);
        assert_eq!(parsed["b"], "hi");
        assert_eq!(body, &payload);
    }

    #[test]
    fn encodes_empty_payload() {
        let m = M { a: 1, b: "".into() };
        let buf = encode_frame(&m, &[]).unwrap();
        assert!(buf.len() >= 5);
        let meta_len = u32::from_le_bytes(buf[..4].try_into().unwrap()) as usize;
        assert_eq!(buf.len(), 4 + meta_len);
    }
}
