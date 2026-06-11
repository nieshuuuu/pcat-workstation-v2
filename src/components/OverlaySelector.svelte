<script lang="ts">
  /**
   * Two-level toggle for selecting which material map to display:
   *   Level 1 — Material (radio chips): CT | Water | Lipid | Total rho
   *   Level 2 — Unit (toggle): Volume % | Mass (mg/mL)
   * Water/lipid only — the 3-material solver (iodine/calcium) was dropped.
   */

  type Props = {
    material: string;
    unit: string;
    onMaterialChange: (material: string) => void;
    onUnitChange: (unit: string) => void;
  };

  let { material, unit, onMaterialChange, onUnitChange }: Props = $props();

  const materials = [
    { key: 'ct', label: 'CT' },
    { key: 'water', label: 'Water' },
    { key: 'lipid', label: 'Lipid' },
    { key: 'density', label: 'Total \u03C1' },
  ];

  const units = [
    { key: 'fraction', label: 'Vol %' },
    { key: 'mass', label: 'mg/mL' },
  ];
</script>

<div class="flex flex-wrap items-center gap-1.5">
  <!-- Material chips -->
  {#each materials as m}
    <button
      class="rounded-full px-2.5 py-1 text-[11px] font-medium transition-colors
             {material === m.key
               ? 'bg-accent text-white'
               : 'bg-surface-tertiary text-text-secondary hover:text-text-primary'}"
      onclick={() => onMaterialChange(m.key)}
    >
      {m.label}
    </button>
  {/each}

  <!-- Divider -->
  <div class="mx-1 h-4 w-px bg-border"></div>

  <!-- Unit toggle pair -->
  {#each units as u}
    <button
      class="rounded-full px-2.5 py-1 text-[11px] font-medium transition-colors
             {unit === u.key
               ? 'bg-accent text-white'
               : 'bg-surface-tertiary text-text-secondary hover:text-text-primary'}
             {material === 'density' || material === 'ct' ? 'opacity-40 pointer-events-none' : ''}"
      onclick={() => onUnitChange(u.key)}
      disabled={material === 'density' || material === 'ct'}
    >
      {u.label}
    </button>
  {/each}
</div>
