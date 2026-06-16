<script lang="ts">
  /**
   * Bottom-right context panel in the 2x2 MPR grid.
   *
   * When pipeline results exist AND centerline is available, shows tabs
   * to switch between CPR view and Analysis dashboard — user never loses
   * access to either.
   */
  import { volumeStore } from '$lib/stores/volumeStore.svelte';
  import CprView from './CprView.svelte';
  import AnalysisDashboard from './AnalysisDashboard.svelte';
  import { pipelineStore } from '$lib/stores/pipelineStore.svelte';
  import { uiStore } from '$lib/stores/uiStore.svelte';

  type Props = {
    phase: 'empty' | 'dicom' | 'seeds' | 'analysis';
  };

  let { phase }: Props = $props();
  let meta = $derived(volumeStore.current);

  // When both CPR and analysis are available, let user switch
  let contextTab = $state<'cpr' | 'analysis'>('cpr');
  let hasResults = $derived(pipelineStore.status === 'complete');
  let hasCpr = $derived(phase === 'seeds' || phase === 'analysis');

  // Fullscreen the Analysis dashboard — it's otherwise confined to this narrow
  // side panel, which is too cramped for the radial/angular charts. The flag
  // lives in uiStore so App's root keydown handler owns Escape: closing the
  // overlay there can't also fire App's "clear vessel" shortcut. (A local
  // window keydown here couldn't prevent that — App's listener is on the same
  // window target and runs regardless of stopPropagation.)
  let maximized = $derived(uiStore.analysisMaximized);

  // Auto-switch to analysis tab when pipeline first completes
  $effect(() => {
    if (pipelineStore.status === 'complete') {
      contextTab = 'analysis';
    }
  });

  // Never leave the overlay stuck open once its results are gone (centerline
  // cleared, patient switched): otherwise the NEXT completed pipeline would pop
  // straight into fullscreen unrequested.
  $effect(() => {
    if (!hasResults) uiStore.analysisMaximized = false;
  });
</script>

<div
  class="flex h-full w-full flex-col bg-surface-secondary text-text-primary"
  class:p-4={phase === 'empty' || phase === 'dicom'}
>
  {#if phase === 'empty'}
    <div class="flex flex-1 flex-col items-center justify-center gap-2">
      <svg
        class="h-10 w-10 text-text-secondary/40"
        fill="none"
        stroke="currentColor"
        stroke-width="1.5"
        viewBox="0 0 24 24"
      >
        <path
          stroke-linecap="round"
          stroke-linejoin="round"
          d="M3.75 9.776c.112-.017.227-.026.344-.026h15.812c.117 0 .232.009.344.026m-16.5 0a2.25 2.25 0 0 0-1.883 2.542l.857 6a2.25 2.25 0 0 0 2.227 1.932H19.05a2.25 2.25 0 0 0 2.227-1.932l.857-6a2.25 2.25 0 0 0-1.883-2.542m-16.5 0V6A2.25 2.25 0 0 1 6 3.75h3.879a1.5 1.5 0 0 1 1.06.44l2.122 2.12a1.5 1.5 0 0 0 1.06.44H18A2.25 2.25 0 0 1 20.25 9v.776"
        />
      </svg>
      <p class="text-sm font-semibold text-text-primary">PCAT Workstation</p>
      <p class="text-xs text-text-secondary">
        Measure pericoronary fat inflammation (FAI) from cardiac CT
      </p>
      <p class="mt-1 text-xs text-text-secondary/60">
        Open a DICOM folder to begin
      </p>
    </div>

  {:else if phase === 'dicom' && meta}
    <div class="flex flex-col gap-3">
      <h3 class="border-b border-border pb-2 text-xs font-semibold tracking-wider text-text-secondary">
        PATIENT INFO
      </h3>
      <div class="flex flex-col gap-2 text-sm">
        <div>
          <span class="text-xs text-text-secondary">Patient Name</span>
          <p class="truncate text-text-primary">{meta.patientName || 'N/A'}</p>
        </div>
        <div>
          <span class="text-xs text-text-secondary">Study Description</span>
          <p class="truncate text-text-primary">{meta.studyDescription || 'N/A'}</p>
        </div>
      </div>
      <h3 class="mt-2 border-b border-border pb-2 text-xs font-semibold tracking-wider text-text-secondary">
        VOLUME
      </h3>
      <div class="flex flex-col gap-2 text-sm">
        <div>
          <span class="text-xs text-text-secondary">Dimensions (Z, Y, X)</span>
          <p class="font-mono text-xs text-text-primary">
            {meta.shape[0]} x {meta.shape[1]} x {meta.shape[2]}
          </p>
        </div>
        <div>
          <span class="text-xs text-text-secondary">Spacing (mm)</span>
          <p class="font-mono text-xs text-text-primary">
            {meta.spacing[0].toFixed(2)} x {meta.spacing[1].toFixed(2)} x {meta.spacing[2].toFixed(2)}
          </p>
        </div>
        <div>
          <span class="text-xs text-text-secondary">Window</span>
          <p class="font-mono text-xs text-text-primary">
            C: {meta.windowCenter} / W: {meta.windowWidth}
          </p>
        </div>
      </div>
    </div>

  {:else if hasCpr}
    <!-- CPR + Analysis: tabbed when results exist -->
    {#if hasResults}
      <!-- Tab bar -->
      <div class="flex shrink-0 border-b border-border bg-surface-secondary">
        <button
          class="px-3 py-1.5 text-xs font-medium transition-colors"
          class:text-accent={contextTab === 'cpr'}
          class:border-b-2={contextTab === 'cpr'}
          class:border-accent={contextTab === 'cpr'}
          class:text-text-secondary={contextTab !== 'cpr'}
          onclick={() => (contextTab = 'cpr')}
        >
          CPR
        </button>
        <button
          class="px-3 py-1.5 text-xs font-medium transition-colors"
          class:text-accent={contextTab === 'analysis'}
          class:border-b-2={contextTab === 'analysis'}
          class:border-accent={contextTab === 'analysis'}
          class:text-text-secondary={contextTab !== 'analysis'}
          onclick={() => (contextTab = 'analysis')}
        >
          Analysis <span class="ml-1 inline-block h-1.5 w-1.5 rounded-full bg-success"></span>
        </button>

        {#if contextTab === 'analysis'}
          <button
            class="ml-auto px-2.5 text-xs font-medium text-text-secondary transition-colors hover:text-accent"
            onclick={() => (uiStore.analysisMaximized = true)}
            title="Fullscreen analysis"
            aria-label="Fullscreen analysis"
          >
            ⤢
          </button>
        {/if}
      </div>
    {/if}

    <!-- Content. When maximized, the dashboard renders in the overlay below
         instead, so it mounts exactly once (one Plotly instance). -->
    <div class="min-h-0 flex-1">
      {#if maximized && hasResults && contextTab === 'analysis'}
        <div class="flex h-full items-center justify-center text-xs text-text-secondary/60">
          Analysis shown fullscreen — press Esc or Close to return.
        </div>
      {:else if hasResults && contextTab === 'analysis'}
        <AnalysisDashboard />
      {:else}
        <CprView />
      {/if}
    </div>
  {/if}
</div>

<!-- Fullscreen Analysis overlay -->
{#if maximized && hasResults}
  <div class="fixed inset-0 z-50 flex flex-col bg-surface">
    <div class="flex shrink-0 items-center justify-between border-b border-border bg-surface-secondary px-4 py-2">
      <span class="text-sm font-semibold text-text-primary">
        FAI Analysis <span class="ml-1 inline-block h-1.5 w-1.5 rounded-full bg-success"></span>
      </span>
      <button
        class="rounded bg-surface-tertiary px-3 py-1 text-xs font-medium text-text-primary hover:bg-surface-tertiary/80"
        onclick={() => (uiStore.analysisMaximized = false)}
        title="Exit fullscreen (Esc)"
      >
        ✕ Close
      </button>
    </div>
    <div class="min-h-0 flex-1">
      <AnalysisDashboard />
    </div>
  </div>
{/if}
