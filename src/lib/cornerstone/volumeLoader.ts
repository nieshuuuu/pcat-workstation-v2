/**
 * Volume loader — constructs a cornerstone3D local volume from a
 * pre-decoded i16 voxel buffer plus the series metadata returned by
 * `loadSeries` in `$lib/api`.
 *
 * No per-slice IPC, no concurrency juggling — the bytes are already here.
 */
import { volumeLoader, cache, type Types } from '@cornerstonejs/core';

import type { VolumeMetadata } from '$lib/api';

type Point3 = Types.Point3;
type Mat3 = Types.Mat3;

/**
 * Build and cache a cornerstone3D volume. Returns the cornerstone volume ID.
 *
 * `volumeKey` is a stable frontend-chosen identifier (e.g., patient+series UID);
 * it is prefixed with `pcat:` to form the cornerstone ID.
 *
 * If a volume with the same ID is already cached, we short-circuit.
 */
export function buildVolume(
  volumeKey: string,
  metadata: VolumeMetadata,
  voxels: Int16Array,
): string {
  const csVolumeId = `pcat:${volumeKey}`;

  const existing = cache.getVolume(csVolumeId);
  if (existing) {
    return csVolumeId;
  }

  // Rust reports [rows, cols, num_slices] (DICOM convention).
  // cornerstone3D uses [X, Y, Z] = [cols, rows, slices].
  const dimensions: Point3 = [metadata.cols, metadata.rows, metadata.num_slices];
  const spacing: Point3 = [
    metadata.pixel_spacing[1],   // sx (column spacing)
    metadata.pixel_spacing[0],   // sy (row spacing)
    metadata.slice_spacing,
  ];
  // Cornerstone takes origin in patient XYZ. Use the full ImagePositionPatient
  // so the world frame matches the DICOM patient frame; the rust side mirrors
  // this in ZYX order. Falls back to slice_positions_z[0] for the Z component
  // when both sides agree (slice_positions_z drives slice sorting).
  const ipp = metadata.image_position_patient;
  const origin: Point3 = [
    ipp[0],
    ipp[1],
    metadata.slice_positions_z[0] ?? ipp[2],
  ];

  // Direction: 3x3 with rows = IOP_row, IOP_col, normal.
  const iopRow: [number, number, number] = [
    metadata.orientation[0], metadata.orientation[1], metadata.orientation[2],
  ];
  const iopCol: [number, number, number] = [
    metadata.orientation[3], metadata.orientation[4], metadata.orientation[5],
  ];
  const normal: [number, number, number] = [
    iopRow[1] * iopCol[2] - iopRow[2] * iopCol[1],
    iopRow[2] * iopCol[0] - iopRow[0] * iopCol[2],
    iopRow[0] * iopCol[1] - iopRow[1] * iopCol[0],
  ];
  const direction = [
    iopRow[0], iopRow[1], iopRow[2],
    iopCol[0], iopCol[1], iopCol[2],
    normal[0], normal[1], normal[2],
  ] as Mat3;

  volumeLoader.createLocalVolume(csVolumeId, {
    scalarData: voxels,
    metadata: {
      BitsAllocated: 16,
      BitsStored: 16,
      SamplesPerPixel: 1,
      HighBit: 15,
      PhotometricInterpretation: 'MONOCHROME2',
      PixelRepresentation: 1,
      Modality: 'CT',
      ImageOrientationPatient: Array.from(direction.slice(0, 6)),
      PixelSpacing: [metadata.pixel_spacing[0], metadata.pixel_spacing[1]],
      FrameOfReferenceUID: `1.2.826.0.1.3680043.8.498.pcat`,
      Columns: metadata.cols,
      Rows: metadata.rows,
      voiLut: [{ windowCenter: metadata.window_center, windowWidth: metadata.window_width }],
      VOILUTFunction: 'LINEAR',
    },
    dimensions,
    spacing,
    origin,
    direction,
  });

  return csVolumeId;
}
