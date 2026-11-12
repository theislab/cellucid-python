"""The viewer's methods, one per action of the catalog: generated, do not edit.

Written by ``scripts/generate_actions.py`` from ``cellucid/action_catalog.json``
(the web app's action catalog, digest 725d32e5357cebbc). Each method sends one
``action`` command through :meth:`AnnDataViewer.act`, which checks and converts
its parameters by the catalog, and returns its :class:`ActionRequest`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING, Any, Literal

from ._actions import OMIT, CellsLike, Omit, ValuesLike

if TYPE_CHECKING:
    from ._viewer import ActionRequest


class ViewerActions:
    """Every action of the catalog as a method (mixed into the viewer)."""

    __slots__ = ()

    if TYPE_CHECKING:

        def act(self, action: str, /, **params: Any) -> ActionRequest: ...

    def analysis_cancel(
        self,
        *,
        analysis: Literal['markers', 'differential', 'dotplot', 'signature', 'composition', 'correlation'],
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Stop an analysis.

        Stop a running analysis, at once, also while its run is still waiting for it; what it was
        showing stays, nothing else changes.

        Action ``analysis.cancel``. It stops work under way in the viewer at once, also while the
        request that started it is still waiting; that work's own step records what it changed, so this
        records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        analysis : {"markers", "differential", "dotplot", "signature", "composition", "correlation"}
            An analysis with a Run button.
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.cancel', analysis=analysis, **{'in': in_})

    def analysis_close(
        self,
    ) -> ActionRequest:
        r"""Close the open analysis.

        Close the analysis open in the Analysis panel (its settings and result stay).

        Action ``analysis.close``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.close')

    def analysis_composition_get(
        self,
        *,
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Read composition options.

        The composition analysis's options, pages and result state.

        Action ``analysis.composition.get``. It reads only and changes nothing; ``request.result()`` is
        what it reads.

        Parameters
        ----------
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.composition.get', **{'in': in_})

    def analysis_composition_set(
        self,
        *,
        options: Mapping[str, Any],
        in_: Literal['panel'] | int = 'panel',
        show: bool = False,
    ) -> ActionRequest:
        r"""Set composition options.

        Change options of the composition analysis (only those given); nothing runs.

        Action ``analysis.composition.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        options : dict
            Keys: ``barsBy`` (optional), ``cellTypes`` (optional), ``cells`` (optional), ``condition``
            (optional), ``correction`` (optional), ``method`` (optional), ``onlyTypes`` (optional),
            ``order`` (optional), ``paired`` (optional), ``samples`` (optional), ``values`` (optional).
        in_ : "panel" or int
            Default ``'panel'``.
        show : bool
            Also open it in the Analysis panel (with in "panel"), as another analysis's "send to" does.
            Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.composition.set', options=options, **{'in': in_}, show=show)

    def analysis_correlation_get(
        self,
        *,
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Read correlation options.

        The correlation analysis's options, pages and result state.

        Action ``analysis.correlation.get``. It reads only and changes nothing; ``request.result()`` is
        what it reads.

        Parameters
        ----------
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.correlation.get', **{'in': in_})

    def analysis_correlation_set(
        self,
        *,
        options: Mapping[str, Any],
        in_: Literal['panel'] | int = 'panel',
        show: bool = False,
    ) -> ActionRequest:
        r"""Set correlation options.

        Change options of the correlation analysis (only those given); nothing runs.

        Action ``analysis.correlation.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        options : dict
            Keys: ``colorBy`` (optional), ``method`` (optional), ``variableX`` (optional), ``variableY``
            (optional).
        in_ : "panel" or int
            Default ``'panel'``.
        show : bool
            Also open it in the Analysis panel (with in "panel"), as another analysis's "send to" does.
            Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.correlation.set', options=options, **{'in': in_}, show=show)

    def analysis_detailed_set(
        self,
        *,
        in_: Literal['panel'] | int = 'panel',
        plot: str | Omit = OMIT,
        plot_options: Sequence[Mapping[str, Any]] | Omit = OMIT,
        variable: None | Mapping[str, Any] | Omit = OMIT,
    ) -> ActionRequest:
        r"""Choose what the detailed analysis plots.

        The variable the detailed analysis plots (a field or a gene, by name, or none), the plot, and
        options of that plot (only those given); it updates by itself.

        Action ``analysis.detailed.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        in_ : "panel" or int
            Default ``'panel'``.
        plot : str
            The plot, by its id as the Plot select offers it for the variable's kind. Optional: left
            out, the viewer decides.
        plot_options : list
            Options of the plot, each by its name and value. Optional: left out, the viewer decides.
        variable : None or dict
            What to plot: a column of the cells (categorical_obs, continuous_obs) or a gene
            (gene_expression) by name, or null for nothing. Keys: ``key``, ``type``. Optional: left out,
            the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.detailed.set', **{'in': in_}, plot=plot, plot_options=plot_options, variable=variable)

    def analysis_differential_get(
        self,
        *,
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Read differential expression options.

        The differential expression analysis's options, pages and result state.

        Action ``analysis.differential.get``. It reads only and changes nothing; ``request.result()`` is
        what it reads.

        Parameters
        ----------
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.differential.get', **{'in': in_})

    def analysis_differential_set(
        self,
        *,
        options: Mapping[str, Any],
        in_: Literal['panel'] | int = 'panel',
        show: bool = False,
    ) -> ActionRequest:
        r"""Set differential expression options.

        Change options of the differential expression analysis (only those given); nothing runs.

        Action ``analysis.differential.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        options : dict
            Keys: ``categoriesA`` (optional), ``categoriesB`` (optional), ``compare`` (optional),
            ``correction`` (optional), ``excludeGenes`` (optional), ``excludePattern`` (optional),
            ``field`` (optional), ``foldChange`` (optional), ``labels`` (optional), ``matrix``
            (optional), ``maxAdjP`` (optional), ``method`` (optional), ``minCells`` (optional),
            ``minCellsPerSample`` (optional), ``minLog2FC`` (optional), ``minPct`` (optional),
            ``onlyGenes`` (optional), ``onlyGenesList`` (optional), ``pConvention`` (optional),
            ``pageA`` (optional), ``pageB`` (optional), ``paired`` (optional), ``replicates``
            (optional), ``within`` (optional), ``withinCategories`` (optional), ``withinField``
            (optional).
        in_ : "panel" or int
            Default ``'panel'``.
        show : bool
            Also open it in the Analysis panel (with in "panel"), as another analysis's "send to" does.
            Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.differential.set', options=options, **{'in': in_}, show=show)

    def analysis_dotplot_get(
        self,
        *,
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Read dot plot options.

        The dot plot analysis's options, pages and result state.

        Action ``analysis.dotplot.get``. It reads only and changes nothing; ``request.result()`` is what
        it reads.

        Parameters
        ----------
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.dotplot.get', **{'in': in_})

    def analysis_dotplot_set(
        self,
        *,
        options: Mapping[str, Any],
        in_: Literal['panel'] | int = 'panel',
        show: bool = False,
    ) -> ActionRequest:
        r"""Set dot plot options.

        Change options of the dot plot analysis (only those given); nothing runs.

        Action ``analysis.dotplot.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        options : dict
            Keys: ``cells`` (optional), ``colorMap`` (optional), ``cutoff`` (optional), ``geneSets``
            (optional), ``groupBy`` (optional), ``groupOrder`` (optional), ``groups`` (optional),
            ``matrix`` (optional), ``maxDotSize`` (optional), ``meanOnlyExpressing`` (optional),
            ``plot`` (optional), ``scale`` (optional), ``swapAxes`` (optional).
        in_ : "panel" or int
            Default ``'panel'``.
        show : bool
            Also open it in the Analysis panel (with in "panel"), as another analysis's "send to" does.
            Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.dotplot.set', options=options, **{'in': in_}, show=show)

    def analysis_export(
        self,
        *,
        analysis: Literal['markers', 'differential', 'dotplot', 'signature', 'composition', 'correlation'],
        all: bool = False,
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Export an analysis result.

        Write an analysis's whole result table as a CSV file, as its expanded view's CSV button does:
        differential expression every gene with scanpy's columns (scores, logfoldchanges, pvals,
        pvals_adj, pct_nz_group, pct_nz_reference, auroc, means, cells); marker genes the shown markers
        of each group, or with all every gene of every group; the dot plot one row per group and gene;
        gene signature every cell's score per set and page (per set and group with a Group by field);
        composition the results per cell type with a test, else the counts; correlation each page's
        statistics. Refused, saying why, when there is no result or it is out of date (run it again
        first).

        Action ``analysis.export``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        analysis : {"markers", "differential", "dotplot", "signature", "composition", "correlation"}
            An analysis with a Run button.
        all : bool
            Marker genes only: every gene of every group (the All genes CSV) instead of the shown
            markers. Default ``False``.
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.export', analysis=analysis, all=all, **{'in': in_})

    def analysis_list(
        self,
    ) -> ActionRequest:
        r"""List the analyses.

        Each analysis of the panel and each analysis window: which is open, its pages, whether it has a
        result, is running, or is stale and why.

        Action ``analysis.list``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.list')

    def analysis_markers_get(
        self,
        *,
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Read marker genes options.

        The marker genes analysis's options, pages and result state.

        Action ``analysis.markers.get``. It reads only and changes nothing; ``request.result()`` is what
        it reads.

        Parameters
        ----------
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.markers.get', **{'in': in_})

    def analysis_markers_set(
        self,
        *,
        options: Mapping[str, Any],
        in_: Literal['panel'] | int = 'panel',
        show: bool = False,
    ) -> ActionRequest:
        r"""Set marker genes options.

        Change options of the marker genes analysis (only those given); nothing runs.

        Action ``analysis.markers.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        options : dict
            Keys: ``cells`` (optional), ``compareAgainst`` (optional), ``correction`` (optional),
            ``direction`` (optional), ``excludeGenes`` (optional), ``excludePattern`` (optional),
            ``foldChange`` (optional), ``groupBy`` (optional), ``matrix`` (optional), ``maxAdjP``
            (optional), ``maxPctRef`` (optional), ``method`` (optional), ``minCells`` (optional),
            ``minLog2FC`` (optional), ``minPctIn`` (optional), ``onlyGenes`` (optional),
            ``onlyGenesList`` (optional), ``pConvention`` (optional), ``rankBy`` (optional),
            ``referenceGroup`` (optional), ``restIncludes`` (optional), ``smallGroups`` (optional),
            ``testGroups`` (optional), ``topN`` (optional).
        in_ : "panel" or int
            Default ``'panel'``.
        show : bool
            Also open it in the Analysis panel (with in "panel"), as another analysis's "send to" does.
            Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.markers.set', options=options, **{'in': in_}, show=show)

    def analysis_open(
        self,
        *,
        analysis: Literal['markers', 'differential', 'dotplot', 'signature', 'composition', 'quick', 'detailed', 'correlation'],
        window: bool = False,
    ) -> ActionRequest:
        r"""Open an analysis.

        Show an analysis in the Analysis panel, or open it in a new floating window with the panel's
        settings.

        Action ``analysis.open``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        analysis : {"markers", "differential", "dotplot", "signature", "composition", "quick", "detailed", "correlation"}
            The analysis: markers (marker genes), differential (differential expression), dotplot,
            signature (gene signature), composition, quick, detailed, correlation.
        window : bool
            Open it in a new floating window. Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.open', analysis=analysis, window=window)

    def analysis_pages_set(
        self,
        *,
        analysis: Literal['markers', 'differential', 'dotplot', 'signature', 'composition', 'quick', 'detailed', 'correlation'],
        pages: Sequence[str],
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Choose the pages of an analysis.

        The highlight pages an analysis compares (by name, "Rest of <page>" for the rest of a page's
        cells); nothing runs.

        Action ``analysis.pages.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        analysis : {"markers", "differential", "dotplot", "signature", "composition", "quick", "detailed", "correlation"}
            The analysis: markers (marker genes), differential (differential expression), dotplot,
            signature (gene signature), composition, quick, detailed, correlation.
        pages : list
            See the action's summary.
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.pages.set', analysis=analysis, pages=pages, **{'in': in_})

    def analysis_quick_set(
        self,
        *,
        categorical: Sequence[str] | Omit = OMIT,
        continuous: Sequence[str] | Omit = OMIT,
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Choose the fields of the quick analysis.

        The categorical columns whose composition and the continuous columns whose statistics the quick
        analysis shows (only the kinds given; before a choice it shows the first two of each); it
        updates by itself.

        Action ``analysis.quick.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        categorical : list
            Categorical columns of the cells, by name (their composition). Optional: left out, the
            viewer decides.
        continuous : list
            Continuous columns of the cells, by name (their statistics). Optional: left out, the viewer
            decides.
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.quick.set', categorical=categorical, continuous=continuous, **{'in': in_})

    def analysis_result(
        self,
        *,
        analysis: Literal['markers', 'differential', 'dotplot', 'signature', 'composition', 'correlation'],
        in_: Literal['panel'] | int = 'panel',
        top: int = 20,
    ) -> ActionRequest:
        r"""Read an analysis result.

        The top rows of an analysis's result table (N a parameter), with how it was tested (groups and
        their cells, test, matrix), what it was computed on and whether it is stale.

        Action ``analysis.result``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        analysis : {"markers", "differential", "dotplot", "signature", "composition", "correlation"}
            An analysis with a Run button.
        in_ : "panel" or int
            Default ``'panel'``.
        top : int
            At most this many rows: per group of marker genes, per direction of a differential
            expression (top 40: 40 up in A and 40 up in B); about.rows says how it counted. Default
            ``20``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.result', analysis=analysis, **{'in': in_}, top=top)

    def analysis_run(
        self,
        *,
        analysis: Literal['markers', 'differential', 'dotplot', 'signature', 'composition', 'correlation'],
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Run an analysis.

        Run an analysis with its options and pages, and wait for its result. The Analysis panel shows
        one analysis at a time: run in it, an analysis is shown first when another is (its section and
        the sidebar must be open).

        Action ``analysis.run``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        analysis : {"markers", "differential", "dotplot", "signature", "composition", "correlation"}
            An analysis with a Run button.
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.run', analysis=analysis, **{'in': in_})

    def analysis_signature_get(
        self,
        *,
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Read gene signature options.

        The gene signature analysis's options, pages and result state.

        Action ``analysis.signature.get``. It reads only and changes nothing; ``request.result()`` is
        what it reads.

        Parameters
        ----------
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.signature.get', **{'in': in_})

    def analysis_signature_set(
        self,
        *,
        options: Mapping[str, Any],
        in_: Literal['panel'] | int = 'panel',
        show: bool = False,
    ) -> ActionRequest:
        r"""Set gene signature options.

        Change options of the gene signature analysis (only those given); nothing runs.

        Action ``analysis.signature.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        options : dict
            Keys: ``bins`` (optional), ``cells`` (optional), ``controlSize`` (optional), ``geneSets``
            (optional), ``groupBy`` (optional), ``groups`` (optional), ``layout`` (optional),
            ``maxRank`` (optional), ``method`` (optional), ``plotType`` (optional), ``rescale``
            (optional), ``seed`` (optional).
        in_ : "panel" or int
            Default ``'panel'``.
        show : bool
            Also open it in the Analysis panel (with in "panel"), as another analysis's "send to" does.
            Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.signature.set', options=options, **{'in': in_}, show=show)

    def analysis_size_set(
        self,
        *,
        analysis: Literal['markers', 'differential', 'dotplot', 'signature', 'composition', 'quick', 'detailed', 'correlation'],
        expanded: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Size an analysis's expanded view.

        The size and splits an analysis's expanded view opens at, remembered per analysis; a person sets
        them by resizing it.

        Action ``analysis.size.set``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        analysis : {"markers", "differential", "dotplot", "signature", "composition", "quick", "detailed", "correlation"}
            The analysis: markers (marker genes), differential (differential expression), dotplot,
            signature (gene signature), composition, quick, detailed, correlation.
        expanded : dict
            The expanded view's box and its splits, in CSS pixels. Keys: ``footerHeight``, ``height``,
            ``optionsWidth``, ``statsWidth``, ``width``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.size.set', analysis=analysis, expanded=expanded)

    def analysis_view_set(
        self,
        *,
        analysis: Literal['markers', 'differential', 'dotplot', 'signature', 'composition', 'correlation'],
        view: Mapping[str, Any],
        in_: Literal['panel'] | int = 'panel',
    ) -> ActionRequest:
        r"""Change how an analysis is shown.

        Open or close its option groups, and choose what its result shows (marker genes: list or
        heatmap, the group listed, the heatmap's look; gene signature: the set shown; the options of the
        plot its expanded view draws); nothing runs.

        Action ``analysis.view.set``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        analysis : {"markers", "differential", "dotplot", "signature", "composition", "correlation"}
            An analysis with an options form.
        view : dict
            What to change of how the analysis is shown (only the keys given). Keys: ``group``
            (optional), ``groups`` (optional), ``heatmap`` (optional), ``more`` (optional), ``plot``
            (optional), ``set`` (optional), ``show`` (optional).
        in_ : "panel" or int
            Default ``'panel'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.view.set', analysis=analysis, view=view, **{'in': in_})

    def analysis_window_close(
        self,
        *,
        window: int,
    ) -> ActionRequest:
        r"""Close an analysis window.

        Close a floating analysis window (its settings and result go with it; undo brings it back).

        Action ``analysis.window.close``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        window : int
            A floating analysis window, by its number (1 is the first open).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.window.close', window=window)

    def analysis_window_collapse(
        self,
        *,
        window: int,
        collapsed: bool,
    ) -> ActionRequest:
        r"""Collapse or expand an analysis window.

        Fold a floating analysis window to its title bar, or unfold it again; its analysis is kept.

        Action ``analysis.window.collapse``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        window : int
            A floating analysis window, by its number (1 is the first open).
        collapsed : bool
            True folds it to its title bar, false unfolds it.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.window.collapse', window=window, collapsed=collapsed)

    def analysis_window_copy(
        self,
        *,
        window: int,
    ) -> ActionRequest:
        r"""Copy an analysis window.

        Open a new floating window with a window's analysis and settings.

        Action ``analysis.window.copy``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        window : int
            A floating analysis window, by its number (1 is the first open).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.window.copy', window=window)

    def analysis_window_place(
        self,
        *,
        window: int,
        left: float,
        top: float,
        width: float,
        height: float,
    ) -> ActionRequest:
        r"""Move or resize an analysis window.

        Give a floating analysis window a place and size; it is kept inside the window.

        Action ``analysis.window.place``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        window : int
            A floating analysis window, by its number (1 is the first open).
        left : float
            See the action's summary.
        top : float
            See the action's summary.
        width : float
            See the action's summary.
        height : float
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('analysis.window.place', window=window, left=left, top=top, width=width, height=height)

    def annotation_publish(
        self,
        *,
        suggestions: Sequence[Mapping[str, Any]] | Omit = OMIT,
        votes: Sequence[Mapping[str, Any]] | Omit = OMIT,
    ) -> ActionRequest:
        r"""Vote and suggest labels.

        Vote suggested labels of categories up or down (none takes the person's vote back) and suggest
        new labels in the person's name (each with their upvote, as the voting dialog's New suggestion
        does), and publish all of them together as one commit: the person's own file in the community
        annotation repository is written as published now plus these changes (a pull request when they
        cannot write to the repository), which others see on their next Pull. Put everything of a turn
        in one call: each call is one commit of the person's file and one question to them. Nothing else
        of the person's is published; votes that change nothing publish nothing. A label already
        suggested for its category is refused with its id (vote on it instead), and nothing of a call is
        published when any of it is refused. Asked of the person first; the history cannot take it back.
        Refused until the person has connected the repository and signed in to GitHub.

        Action ``annotation.publish``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome. It publishes for
        others to see, so the viewer asks the person first, from Python too, as when an agent runs it;
        ``request.result()`` waits for their answer (call it again if its timeout passes first).

        Parameters
        ----------
        suggestions : list
            The new labels, each for one category; one label once per category. Optional: left out, the
            viewer decides.
        votes : list
            The votes, each on one suggestion of one category; one vote per suggestion. Optional: left
            out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('annotation.publish', suggestions=suggestions, votes=votes)

    def annotation_suggestions(
        self,
        *,
        field: Mapping[str, Any],
        category: str | Omit = OMIT,
        cursor: int | Omit = OMIT,
        pull: bool | Omit = OMIT,
    ) -> ActionRequest:
        r"""Read the community annotation.

        The labels suggested for the categories of a field by the community annotation of the dataset
        shown (its connected GitHub repository), with their votes, the person's own vote and each
        category's consensus, as the voting dialog and the legend show them; whether the field is open
        for voting, and when this browser last pulled. The repository is pulled first when this browser
        has not pulled it yet, or when asked (pull), as GitHub sync › Pull latest does: nothing of the
        person's is published. One page holds what fits a tool result; \`next\` starts the following
        one.

        Action ``annotation.suggestions``. It reads only and changes nothing; ``request.result()`` is
        what it reads.

        Parameters
        ----------
        field : dict
            A categorical column of the cells by name (source "obs"). Keys: ``key``, ``source``.
        category : str
            One category; every category of the field when left out. Optional: left out, the viewer
            decides.
        cursor : int
            Where the page starts: the next of the last page; the first page when left out. Optional:
            left out, the viewer decides.
        pull : bool
            Pull the latest suggestions and votes from the repository first (true), or read what this
            browser has (false). Left out: pulled only when this browser has not pulled the repository
            yet. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('annotation.suggestions', field=field, category=category, cursor=cursor, pull=pull)

    def camera_framing(
        self,
        *,
        azimuth: float | Omit = OMIT,
        categories: Sequence[str] | Omit = OMIT,
        elevation: float | Omit = OMIT,
        field: Mapping[str, Any] | Omit = OMIT,
        fill: float = 0.8,
        view: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Frame cells with the camera.

        Where cells sit in a view and cameras that frame them, seen from a direction: the whole view,
        each category of a field (or those named), and those named together. Nothing moves; give a
        camera to camera_set or camera_path_keyframe_add. It frames the view as it is laid out now (side
        by side, a view is narrower): arrange the views first.

        Action ``camera.framing``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        azimuth : float
            Degrees around the vertical axis the camera looks from (3D views); omitted, the view's own.
            90 apart is a quarter turn. Optional: left out, the viewer decides.
        categories : list
            Only these categories of the field; they are also framed together. Omitted, every category
            the view shows. Optional: left out, the viewer decides.
        elevation : float
            Degrees above (positive) or below the subject the camera looks from (3D views); omitted, the
            view's own. Optional: left out, the viewer decides.
        field : dict
            A categorical field whose categories to frame, each on its own. Keys: ``key``, ``source``.
            Optional: left out, the viewer decides.
        fill : float
            How much of the frame the cells fill across or up: 0.8 a close shot with a margin, 0.4 the
            subject with its surroundings. Default ``0.8``.
        view : str
            The view: "View N" as the Compare views list numbers it; the selected view when left out.
            Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera.framing', azimuth=azimuth, categories=categories, elevation=elevation, field=field, fill=fill, view=view)

    def camera_navigation_set(
        self,
        *,
        view: str,
        mode: Literal['orbit', 'planar', 'free'],
    ) -> ActionRequest:
        r"""Set the navigation mode.

        Navigate a view by orbit, planar or free-fly (O, P and G); views with linked cameras share one
        mode.

        Action ``camera.navigation.set``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        view : str
            A view: "View N" as the Compare views list numbers its rows (the live view is View 1 while
            it is shown, the kept views follow in the order they were kept).
        mode : {"orbit", "planar", "free"}
            "orbit", "planar" or "free" (free-fly).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera.navigation.set', view=view, mode=mode)

    def camera_path_clear(
        self,
    ) -> ActionRequest:
        r"""Clear the camera path.

        Delete every keyframe of the camera path; playback stops where the camera is.

        Action ``camera-path.clear``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.clear')

    def camera_path_duration_set(
        self,
        *,
        seconds: float | None,
        keyframe: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Set a camera path transition time.

        Set how long the camera takes from one keyframe to the next, or from every keyframe to its next.

        Action ``camera-path.duration.set``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        seconds : float or None
            Seconds of the transition to the next keyframe, at least 0.1; null paces it automatically.
        keyframe : str
            The keyframe the transition starts at; omitted, every keyframe. Optional: left out, the
            viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.duration.set', seconds=seconds, keyframe=keyframe)

    def camera_path_get(
        self,
    ) -> ActionRequest:
        r"""Read the camera path.

        The camera path: each keyframe (label, camera, when it is reached, how long the camera rests
        there, the time of the transition after it), how it plays, its length, and whether it is
        playing.

        Action ``camera-path.get``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.get')

    def camera_path_hold_set(
        self,
        *,
        seconds: float,
        keyframe: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Set a camera path hold.

        Set how long the camera rests at one keyframe, or at every keyframe, before moving on (0 for no
        rest). The path grows by the hold; each transition keeps its time.

        Action ``camera-path.hold.set``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        seconds : float
            Seconds the camera rests at the keyframe before moving on (a pause in the shot); 0 for none.
        keyframe : str
            The keyframe the camera rests at; omitted, every keyframe. Optional: left out, the viewer
            decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.hold.set', seconds=seconds, keyframe=keyframe)

    def camera_path_keyframe_add(
        self,
        *,
        camera: Mapping[str, Any] | Omit = OMIT,
        embedding: Mapping[str, Any] | None | Omit = OMIT,
        hold: float | Omit = OMIT,
        label: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Add a camera path keyframe.

        Save a camera as the next keyframe of the camera path (before the return to start, when there is
        one), optionally with a hold: seconds the camera rests there.

        Action ``camera-path.keyframe.add``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        camera : dict
            The camera to save, whole (one from camera_framing or views_state); omitted, the selected
            view's camera now. Keys: ``freefly``, ``navigationMode``, ``orbit``. Optional: left out, the
            viewer decides.
        embedding : dict or None
            The embedding the keyframe was recorded on; omitted, the one the selected view shows. Keys:
            ``columns``, ``key``. Optional: left out, the viewer decides.
        hold : float
            Seconds the camera rests at the keyframe before moving on; omitted, 0 (none). Optional: left
            out, the viewer decides.
        label : str
            The keyframe's label, the shot's name: in the Camera path list, and, while the path plays
            with captions on (the default), its caption over the plot from halfway into the move to it
            until halfway out of it, as a documentary's title ("Beta cells", "Overview"; short, a title
            case phrase). Omitted, "KF n", which is never captioned. Unique among the keyframes.
            Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.keyframe.add', camera=camera, embedding=embedding, hold=hold, label=label)

    def camera_path_keyframe_move(
        self,
        *,
        keyframe: str,
        direction: Literal['up', 'down'],
    ) -> ActionRequest:
        r"""Move a camera path keyframe.

        Move a keyframe one place earlier (up) or later (down) in the camera path.

        Action ``camera-path.keyframe.move``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        keyframe : str
            A keyframe, by its label as the Camera path list shows it.
        direction : {"up", "down"}
            "up" one place earlier, "down" one place later.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.keyframe.move', keyframe=keyframe, direction=direction)

    def camera_path_keyframe_remove(
        self,
        *,
        keyframe: str,
    ) -> ActionRequest:
        r"""Delete a camera path keyframe.

        Delete one keyframe of the camera path.

        Action ``camera-path.keyframe.remove``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        keyframe : str
            A keyframe, by its label as the Camera path list shows it.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.keyframe.remove', keyframe=keyframe)

    def camera_path_keyframe_rename(
        self,
        *,
        keyframe: str,
        label: str,
    ) -> ActionRequest:
        r"""Rename a camera path keyframe.

        Give a keyframe another label, unique among the keyframes.

        Action ``camera-path.keyframe.rename``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        keyframe : str
            A keyframe, by its label as the Camera path list shows it.
        label : str
            The new label: the shot's name, and its caption while the path plays with captions on.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.keyframe.rename', keyframe=keyframe, label=label)

    def camera_path_loop_back_set(
        self,
        *,
        enabled: bool,
    ) -> ActionRequest:
        r"""End the camera path at its start.

        Add or remove a last keyframe that returns the camera to the first one (two keyframes or more).

        Action ``camera-path.loop-back.set``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        enabled : bool
            true to end with the return, false to remove it.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.loop-back.set', enabled=enabled)

    def camera_path_play(
        self,
    ) -> ActionRequest:
        r"""Play the camera path.

        Play the camera path in the selected view from its first keyframe, as Play does; answers when it
        ends or camera_path_stop or the person's Stop stops it (the camera stays where it stopped; a
        looping path plays on until stopped, by camera_path_stop or the person's Stop of this turn),
        with ended "end", "stopped", "looping" or "reduced-motion": the person's system asks for reduced
        motion, so the camera went to the last keyframe without moving (seconds: how long it moved, 0
        when at once); nothing played, so never describe a film that ran.

        Action ``camera-path.play``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.play')

    def camera_path_playback_set(
        self,
        *,
        autoplay: bool | Omit = OMIT,
        captions: bool | Omit = OMIT,
        easing: Literal['linear', 'ease-in-out', 'ease-in', 'ease-out'] | Omit = OMIT,
        loop: bool | Omit = OMIT,
        position_interpolation: Literal['catmull-rom', 'linear', 'bezier'] | Omit = OMIT,
        rotation_interpolation: Literal['slerp', 'linear'] | Omit = OMIT,
        speed: int | Omit = OMIT,
    ) -> ActionRequest:
        r"""Set camera path playback.

        Set how the camera path plays: interpolation, easing, speed of automatic transitions, loop,
        autoplay and captions.

        Action ``camera-path.playback.set``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        autoplay : bool
            Start playing as soon as the path has two keyframes. Optional: left out, the viewer decides.
        captions : bool
            Show each named keyframe's label as a caption over the plot while the path plays (on by
            default; a keyframe still called "KF n" has none). Optional: left out, the viewer decides.
        easing : {"linear", "ease-in-out", "ease-in", "ease-out"}
            How the path starts and ends: from rest and to rest (ease-in-out), from rest (ease-in), to
            rest (ease-out), or moving (linear); every transition keeps its own time. Optional: left
            out, the viewer decides.
        loop : bool
            Play the path again from its start when it ends. Optional: left out, the viewer decides.
        position_interpolation : {"catmull-rom", "linear", "bezier"}
            How positions between keyframes are drawn: catmull-rom is the smooth one (no jolt or
            overshoot at a keyframe; a keyframe repeated is a hold). Optional: left out, the viewer
            decides.
        rotation_interpolation : {"slerp", "linear"}
            How the view direction turns between keyframes. Optional: left out, the viewer decides.
        speed : int
            Speed of automatic transitions, 1 (very slow) to 100 (very fast). Optional: left out, the
            viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.playback.set', autoplay=autoplay, captions=captions, easing=easing, loop=loop, position_interpolation=position_interpolation, rotation_interpolation=rotation_interpolation, speed=speed)

    def camera_path_stop(
        self,
    ) -> ActionRequest:
        r"""Stop the camera path.

        Stop a playing camera path where the camera is, at once, also while camera_path_play is still
        waiting for it to end (that play then answers ended "stopped", and its step holds the camera
        where it stopped).

        Action ``camera-path.stop``. It stops work under way in the viewer at once, also while the
        request that started it is still waiting; that work's own step records what it changed, so this
        records no step; ``request.result()`` is its outcome.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera-path.stop')

    def camera_reset(
        self,
        *,
        view: str,
    ) -> ActionRequest:
        r"""Reset the camera.

        Put a view's camera back to its starting place, keeping its navigation mode (the R key).

        Action ``camera.reset``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        view : str
            A view: "View N" as the Compare views list numbers its rows (the live view is View 1 while
            it is shown, the kept views follow in the order they were kept).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera.reset', view=view)

    def camera_set(
        self,
        *,
        view: str,
        camera: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Place the camera.

        Set the camera of a view (views with linked cameras follow); a camera gesture records one of
        these.

        Action ``camera.set``. It changes how the cells are shown, as one step of the viewer's history,
        which the person can undo; ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        view : str
            A view: "View N" as the Compare views list numbers its rows (the live view is View 1 while
            it is shown, the kept views follow in the order they were kept).
        camera : dict
            A camera: its navigation mode, and its orbit and free-fly poses (the one the mode drives is
            used). Coordinates are the view's (the embedding scaled to about -1..1, y up); an orbit
            camera sits at target + radius·(sin phi·cos theta, cos phi, sin phi·sin theta) looking at
            target (theta the azimuth, phi from straight above, radians); a planar one at target + (0,
            0, radius). Take cameras from \`camera_framing\` rather than composing them. Keys:
            ``freefly``, ``navigationMode``, ``orbit``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('camera.set', view=view, camera=camera)

    def category_delete(
        self,
        *,
        field: Mapping[str, Any],
        category: str,
        unassigned_label: str = 'Unassigned',
    ) -> ActionRequest:
        r"""Delete a category.

        Move every cell of one category into the Unassigned category. A column from the file becomes a
        new column (the original moves to Deleted fields); a column made in the app is edited in place.

        Action ``category.delete``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A categorical column of the cells by name (source "obs"). Keys: ``key``, ``source``.
        category : str
            A category label as the legend shows it.
        unassigned_label : str
            The category the cells move into, made when missing. Default ``'Unassigned'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('category.delete', field=field, category=category, unassigned_label=unassigned_label)

    def category_merge(
        self,
        *,
        field: Mapping[str, Any],
        category: str,
        into: str,
        label: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Merge two categories.

        Move every cell of one category into another. The merged category is named "merged <category> +
        <into>" (as the legend's drag names it) unless label names it: give label the into category's
        own name to keep it. A column from the file becomes a new column (the original moves to Deleted
        fields); a column made in the app is edited in place.

        Action ``category.merge``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A categorical column of the cells by name (source "obs"). Keys: ``key``, ``source``.
        category : str
            The category whose cells move.
        into : str
            The category they join.
        label : str
            The merged category's name; "merged <category> + <into>" when left out. Optional: left out,
            the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('category.merge', field=field, category=category, into=into, label=label)

    def category_rename(
        self,
        *,
        field: Mapping[str, Any],
        category: str,
        name: str,
    ) -> ActionRequest:
        r"""Rename a category.

        Give one category of a categorical column another label; its original label stays known.

        Action ``category.rename``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A categorical column of the cells by name (source "obs"). Keys: ``key``, ``source``.
        category : str
            A category label as the legend shows it.
        name : str
            The new label, unique in the column.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('category.rename', field=field, category=category, name=name)

    def cells_export(
        self,
        *,
        format: Literal['csv', 'h5ad'] = 'csv',
        page: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Export the highlighted cells.

        Write the highlighted cells of a highlight page as a file, as the Highlights list's export
        buttons do: an obs-column CSV (every open cell, a true/false column per highlight, which
        cellucid.pull() in Python and cellucid_pull() in R read back) or an .h5ad subset of the
        highlighted cells with their expression rows. The .h5ad reads the cells' expression first; its
        time and bytes are estimated from this connection before it starts.

        Action ``cells.export``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        format : {"csv", "h5ad"}
            csv: the obs-column CSV; h5ad: an AnnData file of the highlighted cells. Default ``'csv'``.
        page : str
            The highlight page; the page shown when left out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('cells.export', format=format, page=page)

    def cells_open(
        self,
        *,
        mode: Literal['all', 'subsample'],
        allocation: Literal['proportional', 'balanced'] = 'proportional',
        count: int | Omit = OMIT,
        embedding: str | Omit = OMIT,
        field: str | Omit = OMIT,
        fraction: float | Omit = OMIT,
        method: Literal['uniform', 'blocks', 'stratified', 'spatial'] | None = None,
        seed: int | None = None,
    ) -> ActionRequest:
        r"""Open a share of the cells.

        Read the dataset again with every cell, or with a subsample: an exact \`count\` of cells or a
        \`fraction\` of them, chosen by a subsample method; highlights, labels, views and filters are
        kept for the cells still open.

        Action ``cells.open``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        mode : {"all", "subsample"}
            all: every cell; subsample: the count or fraction given.
        allocation : {"proportional", "balanced"}
            Stratified: how the cells are split over the categories. Default ``'proportional'``.
        count : int
            How many cells, exactly (at least 2; the dataset's count or more opens every cell). Not with
            fraction. Optional: left out, the viewer decides.
        embedding : str
            Spatial: the obsm embedding it thins. Optional: left out, the viewer decides.
        field : str
            Stratified: the categorical column it is stratified by. Optional: left out, the viewer
            decides.
        fraction : float
            The share of the dataset's cells: 0.1 opens a tenth, 1 every cell. Not with count. Optional:
            left out, the viewer decides.
        method : {"uniform", "blocks", "stratified", "spatial", null}
            How the cells are chosen; null for the file's default. Default ``None``.
        seed : int or None
            The subsample's seed; null for the file's. Default ``None``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('cells.open', mode=mode, allocation=allocation, count=count, embedding=embedding, field=field, fraction=fraction, method=method, seed=seed)

    def cells_read_genes(
        self,
        *,
        cells: int,
    ) -> ActionRequest:
        r"""Read genes for more cells.

        Read every gene for the first \`cells\` open cells (a dataset stored cell by cell reads genes
        for a sample of the open cells); stopping keeps what was read.

        Action ``cells.read-genes``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        cells : int
            How many of the open cells; more than are open reads every one.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('cells.read-genes', cells=cells)

    def cells_values(
        self,
        *,
        cells: CellsLike,
        fields: Sequence[Mapping[str, Any]],
    ) -> ActionRequest:
        r"""Read values of cells.

        The values of named cells (original positions) in each field or gene asked for.

        Action ``cells.values``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        cells : cells (rows, a boolean mask or obs names)
            Cells by original position (their row in the file's obs).
        fields : list
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('cells.values', cells=cells, fields=fields)

    def cells_visibility_set(
        self,
        *,
        cells: CellsLike | None,
        visible: bool,
    ) -> ActionRequest:
        r"""Hide or show cells.

        Hide or show cells one by one in the selected view, on top of every filter.

        Action ``cells.visibility.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        cells : cells (rows, a boolean mask or obs names)
            Cells by original position, or null for every cell.
        visible : bool
            true to show them, false to hide them.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('cells.visibility.set', cells=cells, visible=visible)

    def color_by(
        self,
        *,
        field: Mapping[str, Any] | None,
    ) -> ActionRequest:
        r"""Color the cells by a field or gene.

        Color the view in focus by a column of the cells (source "obs") or a gene (source "var"),
        reading it first when needed; null clears the coloring.

        Action ``color.by``. It changes how the cells are shown, as one step of the viewer's history,
        which the person can undo; ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict or None
            The field or gene by name; null colors the cells by nothing. Keys: ``key``, ``source``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('color.by', field=field)

    def color_category(
        self,
        *,
        field: Mapping[str, Any],
        category: str,
        color: str,
    ) -> ActionRequest:
        r"""Color one category.

        Give one category of a categorical column another color in the view in focus: each view keeps
        its own colors, so view_focus first selects the view to change (a view not colored by the column
        shows no change).

        Action ``color.category``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        field : dict
            A categorical column of the cells by name (source "obs"). Keys: ``key``, ``source``.
        category : str
            A category label as the legend shows it.
        color : str
            A color as six hex digits, "#e41a1c".

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('color.category', field=field, category=category, color=color)

    def color_continuous_colormap(
        self,
        *,
        field: Mapping[str, Any],
        colormap: Literal['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'turbo', 'twilight', 'coolwarm', 'bluered', 'spectral', 'rdylbu', 'rdbu', 'piyg', 'prgn', 'brbg', 'seismic', 'ylgnbu', 'ylorrd', 'ylorbr', 'gnbu', 'bupu', 'rdpu', 'purd', 'orrd', 'blues', 'greens', 'reds', 'purples', 'oranges', 'greys', 'heat', 'hot', 'cool', 'jet', 'rainbow', 'ocean', 'terrain', 'spring', 'summer', 'autumn', 'winter', 'bone', 'copper', 'pink'],
    ) -> ActionRequest:
        r"""Choose a color map.

        The color map a continuous field or gene is drawn with in the view in focus (each view keeps its
        own).

        Action ``color.continuous.colormap``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The continuous field or gene the view in focus is colored by, by name. Keys: ``key``,
            ``source``.
        colormap : {"viridis", "plasma", "inferno", "magma", "cividis", "turbo", "twilight", "coolwarm", "bluered", "spectral", "rdylbu", "rdbu", "piyg", "prgn", "brbg", "seismic", "ylgnbu", "ylorrd", "ylorbr", "gnbu", "bupu", "rdpu", "purd", "orrd", "blues", "greens", "reds", "purples", "oranges", "greys", "heat", "hot", "cool", "jet", "rainbow", "ocean", "terrain", "spring", "summer", "autumn", "winter", "bone", "copper", "pink"}
            The color map by its id (viridis, magma, coolwarm, …).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('color.continuous.colormap', field=field, colormap=colormap)

    def color_continuous_log(
        self,
        *,
        field: Mapping[str, Any],
        enabled: bool,
    ) -> ActionRequest:
        r"""Log color scale.

        Color a continuous field or gene on a log scale, or a linear one, in the view in focus (each
        view keeps its own).

        Action ``color.continuous.log``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        field : dict
            The continuous field or gene the view in focus is colored by, by name. Keys: ``key``,
            ``source``.
        enabled : bool
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('color.continuous.log', field=field, enabled=enabled)

    def color_continuous_range(
        self,
        *,
        field: Mapping[str, Any],
        min: float,
        max: float,
    ) -> ActionRequest:
        r"""Set the color range.

        The values the color map spans for a continuous field or gene that is rescaled to a range, in
        the view in focus (each view keeps its own).

        Action ``color.continuous.range``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The continuous field or gene the view in focus is colored by, by name. Keys: ``key``,
            ``source``.
        min : float
            See the action's summary.
        max : float
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('color.continuous.range', field=field, min=min, max=max)

    def color_continuous_rescale(
        self,
        *,
        field: Mapping[str, Any],
        enabled: bool,
    ) -> ActionRequest:
        r"""Rescale colors to the visible range.

        Stretch the color map over the visible range of a continuous field or gene, or over all its
        values, in the view in focus (each view keeps its own).

        Action ``color.continuous.rescale``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The continuous field or gene the view in focus is colored by, by name. Keys: ``key``,
            ``source``.
        enabled : bool
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('color.continuous.rescale', field=field, enabled=enabled)

    def color_continuous_zero_gray(
        self,
        *,
        field: Mapping[str, Any],
        enabled: bool,
    ) -> ActionRequest:
        r"""Gray for zero.

        Draw a gene's zeros in gray beneath the expressing cells, or in the color map's lowest color, in
        the view in focus (each view keeps its own).

        Action ``color.continuous.zero-gray``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The continuous field or gene the view in focus is colored by, by name. Keys: ``key``,
            ``source``.
        enabled : bool
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('color.continuous.zero-gray', field=field, enabled=enabled)

    def dataset_open(
        self,
        *,
        name: str,
        allocation: Literal['proportional', 'balanced'] = 'proportional',
        count: int | Omit = OMIT,
        embedding: str | Omit = OMIT,
        field: str | Omit = OMIT,
        fraction: float | Omit = OMIT,
        kind: Literal['any', 'sample', 'my-datasets', 'link'] = 'any',
        method: Literal['uniform', 'blocks', 'stratified', 'spatial'] | None = None,
        seed: int | None = None,
    ) -> ActionRequest:
        r"""Open a dataset.

        Open a dataset by the name the person gives: a built-in sample, one of their datasets (My
        datasets), a link, or an accession the link box takes. It replaces the open dataset; that one's
        history stays kept in this browser while Autosave keeps it (steps it does not keep are asked
        about first) and the new one starts its own. A link, accession or server the person has not
        typed, pasted or opened in this tab is asked of them first (before anything is requested), and
        so is a whole-file download. A count or fraction opens that subsample of its cells at once, read
        once (as cells_open takes it); without, it opens with the cells the person chose for openings
        (Session › Subsampling), or the cells one of My datasets was last open with.

        Action ``dataset.open``. It opens another dataset in place of the one the viewer shows, as an
        agent in the viewer does. A notebook's viewer shows the dataset its notebook serves, so the
        viewer refuses it from Python: show another dataset with :func:`cellucid.show_anndata` or
        :func:`cellucid.serve_anndata`.

        Parameters
        ----------
        name : str
            A sample's name or id, a name from My datasets, a link to a dataset (a file, a Zarr store, a
            dataset folder, a CELLxGENE, Zenodo, Figshare or GitHub address), an accession (GSE123456, a
            DOI), or a Cellucid server (localhost:8765).
        allocation : {"proportional", "balanced"}
            Stratified: how the cells are split over the categories. Default ``'proportional'``.
        count : int
            How many cells, exactly (at least 2; the dataset's count or more opens every cell). Not with
            fraction. Optional: left out, the viewer decides.
        embedding : str
            Spatial: the obsm embedding it thins. Optional: left out, the viewer decides.
        field : str
            Stratified: the categorical column it is stratified by. Optional: left out, the viewer
            decides.
        fraction : float
            The share of the dataset's cells: 0.1 opens a tenth, 1 every cell. Not with count. Optional:
            left out, the viewer decides.
        kind : {"any", "sample", "my-datasets", "link"}
            Where to look for the name: any (samples, then My datasets, then as a link), sample,
            my-datasets or link. Default ``'any'``.
        method : {"uniform", "blocks", "stratified", "spatial", null}
            How the cells are chosen; null for the file's default. Default ``None``.
        seed : int or None
            The subsample's seed; null for the file's. Default ``None``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('dataset.open', name=name, allocation=allocation, count=count, embedding=embedding, field=field, fraction=fraction, kind=kind, method=method, seed=seed)

    def dataset_share_link(
        self,
    ) -> ActionRequest:
        r"""Make a share link.

        The public link that opens the dataset shown on www.cellucid.com for anyone, with its "Open in
        Cellucid" badge and an embed snippet, as Share in Session › Open data makes them; what a link
        carried to grant access is left out. Refused for a dataset no link reopens (a file of this
        computer, a notebook's, a server on this computer).

        Action ``dataset.share-link``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome. It publishes for
        others to see, so the viewer asks the person first, from Python too, as when an agent runs it;
        ``request.result()`` waits for their answer (call it again if its timeout passes first).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('dataset.share-link')

    def dataset_summary(
        self,
    ) -> ActionRequest:
        r"""Read the dataset summary.

        The open dataset: its name, description, citation, how many cells (and how many are open),
        genes, embeddings, vector fields and fields.

        Action ``dataset.summary``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('dataset.summary')

    def datasets_list(
        self,
    ) -> ActionRequest:
        r"""List the datasets to open.

        The datasets this tab can open by name: the built-in samples and the person's own datasets
        opened before in this browser (My datasets, each with its size and what this browser keeps of
        it), with the one open now.

        Action ``datasets.list``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('datasets.list')

    def display_reset(
        self,
    ) -> ActionRequest:
        r"""Reset the camera and display.

        The Reset camera button: the selected view's camera, its dimension's navigation mode, the
        background, the render mode and the smoke resolutions back to their defaults.

        Action ``display.reset``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('display.reset')

    def expression_transform_set(
        self,
        *,
        log1p: bool | Omit = OMIT,
        normalize: bool | Omit = OMIT,
    ) -> ActionRequest:
        r"""Set the expression values.

        Show and analyze the genes normalized to 10,000 counts per cell and/or log1p (Data mapping ›
        Expression values; scanpy normalize_total(target_sum=1e4) then log1p). One setting for coloring,
        legends, filters, the cell card, reads and every analysis (pseudobulk differential expression
        always sums raw counts). The person's setting: an agent changes it only when the person asks.
        Offered when the mapped matrix holds counts; normalizing reads each cell's total from an obs
        column, else sums every gene once (heavy on a large matrix: asked first).

        Action ``expression.transform.set``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        log1p : bool
            ln(1 + x), after the normalization when both are on. Optional: left out, the viewer decides.
        normalize : bool
            Each cell's values scaled to 10,000 counts: x · 10,000 / the cell's total. Optional: left
            out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('expression.transform.set', log1p=log1p, normalize=normalize)

    def field_categories(
        self,
        *,
        field: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Read the categories of a field.

        Each category of a categorical field with its cells, the cells shown in the active view, its
        color and whether it is shown.

        Action ``field.categories``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        field : dict
            A field by name: source "obs" for a column of the cells, "var" for a gene; key is the name
            as shown. Keys: ``key``, ``source``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.categories', field=field)

    def field_consensus(
        self,
        *,
        key: str,
        field: Mapping[str, Any],
        mapping: Sequence[Mapping[str, Any]],
        settings: Mapping[str, Any],
        provenance: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Build the consensus column.

        The community annotation's consensus for a categorical column as a column of its own: each cell
        gets the agreed label of its category, or the status of a category without one (disputed,
        pending); built again under the same name, it is replaced. The outcome of each category is named
        in the parameters, so it replays by name on any dataset with the same column and categories.

        Action ``field.consensus``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        key : str
            The consensus column's name.
        field : dict
            The annotated column the consensus is of. Keys: ``key``, ``source``.
        mapping : list
            The outcome of each category of the column (a category left out is pending).
        settings : dict
            The consensus rule it was built with. Keys: ``minAnnotators``, ``threshold``.
        provenance : dict
            Where the votes came from. Keys: ``branch``, ``commit``, ``pulledAt``, ``repository``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.consensus', key=key, field=field, mapping=mapping, settings=settings, provenance=provenance)

    def field_create_from_pages(
        self,
        *,
        key: str,
        pages: Sequence[Mapping[str, Any]],
        intersection_labels: Sequence[Mapping[str, Any]] | Omit = OMIT,
        overlap: Literal['first', 'last', 'overlap-label', 'intersections'] = 'first',
        overlap_label: str = 'Overlap',
        rest: Mapping[str, Any] | Omit = OMIT,
        uncovered_label: str = 'Unassigned',
    ) -> ActionRequest:
        r"""Make a field from highlight pages.

        A new categorical column whose categories are highlight pages, shown at once: each page's cells
        take its label; cells on several pages follow the overlap rule; cells on none keep their
        category of rest when it is given (refine an annotation: a page per new label, every other cell
        keeps its label), else take uncovered_label, a category only when there are such cells.

        Action ``field.create-from-pages``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        key : str
            The new column's name.
        pages : list
            The pages, in order (first and last refer to it), each with its category label.
        intersection_labels : list
            With overlap "intersections": the label of each combination of pages that has cells.
            Optional: left out, the viewer decides.
        overlap : {"first", "last", "overlap-label", "intersections"}
            Cells on several pages: the first page's label, the last's, the overlap label, or one label
            per combination of pages. Default ``'first'``.
        overlap_label : str
            The label of cells on several pages when overlap is "overlap-label". Default ``'Overlap'``.
        rest : dict
            A categorical column whose category the cells on no page keep (split one label: a page per
            new label, every other cell keeps its label in one call). A page label that is one of its
            categories is that category. Keys: ``key``, ``source``. Optional: left out, the viewer
            decides.
        uncovered_label : str
            The category of cells on no page (with rest: of those with no category there); empty leaves
            them without one. Default ``'Unassigned'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.create-from-pages', key=key, pages=pages, intersection_labels=intersection_labels, overlap=overlap, overlap_label=overlap_label, rest=rest, uncovered_label=uncovered_label)

    def field_delete(
        self,
        *,
        field: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Delete a field.

        Remove a column of the cells or a gene from the lists; it can be restored from Deleted fields
        until purged.

        Action ``field.delete``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A field by name: source "obs" for a column of the cells, "var" for a gene; key is the name
            as shown. Keys: ``key``, ``source``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.delete', field=field)

    def field_duplicate(
        self,
        *,
        field: Mapping[str, Any],
        name: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Duplicate a field.

        Copy a column of the cells or a gene as a new field (named "<name> (copy)" unless a name is
        given) and color by the copy.

        Action ``field.duplicate``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A field by name: source "obs" for a column of the cells, "var" for a gene; key is the name
            as shown. Keys: ``key``, ``source``.
        name : str
            The copy's name; by default "<name> (copy)", numbered when taken. Optional: left out, the
            viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.duplicate', field=field, name=name)

    def field_export(
        self,
        *,
        field: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Export a field.

        Write a column of the cells (categorical or continuous) as an obs-column CSV, as its Export
        button in Color by does: every open cell and its value, which cellucid.pull() in Python and
        cellucid_pull() in R read back.

        Action ``field.export``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        field : dict
            A column of the cells by name (source "obs"), categorical or continuous. Keys: ``key``,
            ``source``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.export', field=field)

    def field_group(
        self,
        *,
        field: Mapping[str, Any],
        groups: Sequence[Mapping[str, Any]],
        name: str | Omit = OMIT,
        rest: str | None = None,
    ) -> ActionRequest:
        r"""Group categories into a new field.

        A new categorical column from another, whose categories are groups of its categories (a lineage
        from cell types: Ductal and Ngn3 low EP as "Progenitor", Alpha, Beta, Delta and Epsilon as
        "Endocrine"), shown at once; the column it is made from stays. A category no group lists is kept
        as it is, or joins rest when rest names a category. One step: one undo.

        Action ``field.group``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The categorical column whose categories are grouped. Keys: ``key``, ``source``.
        groups : list
            The new categories, in order, each with the categories it holds (a category in one group at
            most).
        name : str
            The new column's name; by default "<name> (grouped)", numbered when taken. Optional: left
            out, the viewer decides.
        rest : str or None
            The category every category no group lists joins ("Other", or a group's label); null keeps
            each as it is. Default ``None``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.group', field=field, groups=groups, name=name, rest=rest)

    def field_purge(
        self,
        *,
        field: Mapping[str, Any],
        which: int | Omit = OMIT,
    ) -> ActionRequest:
        r"""Delete a deleted field for good.

        Delete a deleted field for good (Deleted fields › Delete for good…): it leaves Deleted fields
        and can no longer be restored there (undo still brings it back).

        Action ``field.purge``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The deleted field by the name it had when deleted. Keys: ``key``, ``source``.
        which : int
            Only when several deleted fields have this name: which, 1 the first of them as Deleted
            fields lists them. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.purge', field=field, which=which)

    def field_push(
        self,
        *,
        key: str,
        kind: Literal['categorical', 'continuous'],
        categories: Sequence[str] | None,
        n_obs: int,
        cells: CellsLike | None,
        values: ValuesLike,
    ) -> ActionRequest:
        r"""Add a column of values.

        A notebook's column as a field of the cells: categorical (codes into its categories) or
        continuous; a later push of the same name replaces it.

        Action ``field.push``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        key : str
            The column's name.
        kind : {"categorical", "continuous"}
            See the action's summary.
        categories : list or None
            A categorical column's categories (its values are codes into them; the largest code is "no
            category"); null for a continuous one.
        n_obs : int
            The cells of the dataset the column belongs to.
        cells : cells (rows, a boolean mask or obs names)
            The cells the values are of, as original positions; null for every cell in order.
        values : array of numbers, one per cell
            One value per cell named: a code (categorical) or a number (continuous).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.push', key=key, kind=kind, categories=categories, n_obs=n_obs, cells=cells, values=values)

    def field_rename(
        self,
        *,
        field: Mapping[str, Any],
        name: str,
    ) -> ActionRequest:
        r"""Rename a field.

        Give a column of the cells or a gene another name; its original name stays known, and renaming
        it back to that name restores it.

        Action ``field.rename``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A field by name: source "obs" for a column of the cells, "var" for a gene; key is the name
            as shown. Keys: ``key``, ``source``.
        name : str
            A field name: not empty, unique among the fields shown.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.rename', field=field, name=name)

    def field_restore(
        self,
        *,
        field: Mapping[str, Any],
        which: int | Omit = OMIT,
    ) -> ActionRequest:
        r"""Restore a deleted field.

        Bring back a deleted field under its name, as it was when deleted.

        Action ``field.restore``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The deleted field by the name it had when deleted. Keys: ``key``, ``source``.
        which : int
            Only when several deleted fields have this name: which, 1 the first of them as Deleted
            fields lists them. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.restore', field=field, which=which)

    def field_revert_name(
        self,
        *,
        field: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Restore a field's original name.

        Give a renamed field back the name it has in the file.

        Action ``field.revert-name``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A field by name: source "obs" for a column of the cells, "var" for a gene; key is the name
            as shown. Keys: ``key``, ``source``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.revert-name', field=field)

    def field_summary(
        self,
        *,
        field: Mapping[str, Any],
        by: Mapping[str, Any],
        categories: Sequence[str] | Omit = OMIT,
        top: int = 50,
        visible_only: bool = False,
        within: Mapping[str, Any] | Omit = OMIT,
    ) -> ActionRequest:
        r"""Summarize a field per category.

        A continuous field or a gene per category of a categorical field (cells, mean, median, range,
        share above zero), or a categorical field's top categories per category.

        Action ``field.summary``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        field : dict
            A field by name: source "obs" for a column of the cells, "var" for a gene; key is the name
            as shown. Keys: ``key``, ``source``.
        by : dict
            The categorical obs field to group the cells by. Keys: ``key``, ``source``.
        categories : list
            Only these categories of \`by\`, by name (else every category with cells). Optional: left
            out, the viewer decides.
        top : int
            For a categorical field: at most this many of its categories per group. Default ``50``.
        visible_only : bool
            Only the cells the active view shows. Default ``False``.
        within : dict
            Only the cells of some categories of another categorical field, or in a range of a
            continuous field or gene (both ends included), leaving the view as it is: donors within one
            cell type, cell types per age range. Keys: ``categories`` (optional), ``field``, ``max``
            (optional), ``min`` (optional). Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('field.summary', field=field, by=by, categories=categories, top=top, visible_only=visible_only, within=within)

    def fields_list(
        self,
        *,
        kind: Literal['categorical', 'continuous'] | Omit = OMIT,
        query: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""List the fields.

        The cell annotations (obs fields) with their kind and, for a categorical one, how many
        categories it has.

        Action ``fields.list``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        kind : {"categorical", "continuous"}
            Only fields of this kind. Optional: left out, the viewer decides.
        query : str
            Only fields whose key contains this text (any case). Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('fields.list', kind=kind, query=query)

    def figure_export(
        self,
    ) -> ActionRequest:
        r"""Export the figure.

        Write the figure of the view shown (every view in a grid when the settings say so) as its
        Download button does, with the figure export's settings as they are (figure_export_set changes
        them: format, size, labels, legend, style): one SVG or PNG file, or several in one ZIP archive
        (view.picture with source "export" shows that figure first, nothing written). Refused, saying
        why, when the figure cannot reproduce what the view shows (a render mode other than points, an
        overlay it cannot draw).

        Action ``figure.export``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('figure.export')

    def figure_export_get(
        self,
    ) -> ActionRequest:
        r"""Read the figure export settings.

        The Figure export panel's settings now, each named as figure-export.set names it: what Export
        would write, and what a session or the history put back.

        Action ``figure-export.get``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('figure-export.get')

    def figure_export_set(
        self,
        *,
        auto_text_size: bool | Omit = OMIT,
        background: Literal['viewer', 'white', 'transparent', 'custom'] | Omit = OMIT,
        background_color: str | Omit = OMIT,
        depth_sort: bool | Omit = OMIT,
        download: Literal['svg', 'png', 'svg+png', 'png-multi', 'all'] | Omit = OMIT,
        dpi: Literal[150, 300, 600] | Omit = OMIT,
        emphasize_selection: bool | Omit = OMIT,
        export_all_views: bool | Omit = OMIT,
        font_family: Literal['Arial, Helvetica, sans-serif', 'Helvetica, Arial, sans-serif', '"Times New Roman", Times, serif'] | Omit = OMIT,
        font_sizes: Mapping[str, Any] | None | Omit = OMIT,
        frame: Mapping[str, Any] | None | Literal['cells'] | Omit = OMIT,
        include_axes: bool | Omit = OMIT,
        include_legend: bool | Omit = OMIT,
        legend_position: Literal['right', 'bottom'] | Omit = OMIT,
        optimized_target_count: int | Omit = OMIT,
        preview: Mapping[str, Any] | Omit = OMIT,
        reference_grid: bool | Omit = OMIT,
        selection_muted_opacity: float | Omit = OMIT,
        show_orientation: bool | Omit = OMIT,
        size_preset: Literal['single-column', 'one-and-half-column', 'double-column', 'custom'] | Omit = OMIT,
        svg_strategy: Literal['full-vector', 'optimized-vector', 'hybrid'] | Omit = OMIT,
        title: str | None | Omit = OMIT,
        title_follows_dataset: bool | Omit = OMIT,
        width: float | Omit = OMIT,
        x_label: str | Omit = OMIT,
        y_label: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Change figure export settings.

        Set Figure export options: width (millimetres or a journal column; the figure fills it, text in
        points keeps its printed size), frame (frame "cells" fits it to the cells, so plots take their
        shape, not the pane's), background, fonts, title and axis labels, axes, legend, preview, what
        Export downloads; only those given change, nothing is exported (figure-export.get reads them;
        view.picture with source "export" shows the figure they make).

        Action ``figure-export.set``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        auto_text_size : bool
            Set the text in the sizes journals print figure text at: 7 pt, ticks and centroid labels 6
            pt, the title 8 pt (font_sizes is then null). Optional: left out, the viewer decides.
        background : {"viewer", "white", "transparent", "custom"}
            The background: the viewer's, white, transparent or custom. Optional: left out, the viewer
            decides.
        background_color : str
            The custom background color, lowercase #rrggbb. Optional: left out, the viewer decides.
        depth_sort : bool
            Draw 3D points back to front. Optional: left out, the viewer decides.
        download : {"svg", "png", "svg+png", "png-multi", "all"}
            What Export downloads. Optional: left out, the viewer decides.
        dpi : {150, 300, 600}
            The resolution of a PNG: its pixels are the width in inches times this (an SVG is the same
            size, as vectors). Optional: left out, the viewer decides.
        emphasize_selection : bool
            Fade the cells that are not highlighted. Optional: left out, the viewer decides.
        export_all_views : bool
            Export every view of the multiview, side by side: each panel is titled by its view's name
            (view_rename; a view with none by its field and filters). Optional: left out, the viewer
            decides.
        font_family : {"Arial, Helvetica, sans-serif", "Helvetica, Arial, sans-serif", "\\"Times New Roman\\", Times, serif"}
            The font of every text. Optional: left out, the viewer decides.
        font_sizes : dict or None
            Text sizes in points (1/72 inch; halves allowed) while auto_text_size is off, else null.
            They keep their printed size at any width and resolution. Keys: ``axis_label``, ``base``,
            ``centroid_label``, ``legend``, ``ticks``, ``title``. Optional: left out, the viewer
            decides.
        frame : dict or None or "cells"
            The part of the view exported (needs the preview on, which "cells" turns on), null for the
            whole view, or "cells": the frame every cell the figure draws fills, with a small margin
            (the panel's Fit the cells), so a plot takes the cells' shape instead of its pane's (a
            narrow pane gives a narrow plot; never rearrange the person's screen for a figure). With
            export_all_views, one frame holding every panel's cells. Keys: ``lock_aspect``, ``rect``.
            Optional: left out, the viewer decides.
        include_axes : bool
            Draw the axes. Optional: left out, the viewer decides.
        include_legend : bool
            Draw the legend. Optional: left out, the viewer decides.
        legend_position : {"right", "bottom"}
            Where the legend goes. Optional: left out, the viewer decides.
        optimized_target_count : int
            How many points an optimized SVG keeps. Optional: left out, the viewer decides.
        preview : dict
            The preview in the panel and the color-vision deficiency it simulates; with the preview on,
            view.picture with source "export" draws that simulation too (Export writes the true colors).
            Keys: ``colorblind_mode``, ``enabled``. Optional: left out, the viewer decides.
        reference_grid : bool
            Draw a reference grid. Optional: left out, the viewer decides.
        selection_muted_opacity : float
            The opacity of the faded cells, 0 to 1. Optional: left out, the viewer decides.
        show_orientation : bool
            Draw the 3D orientation indicator. Optional: left out, the viewer decides.
        size_preset : {"single-column", "one-and-half-column", "double-column", "custom"}
            The figure's width as a journal column: single-column (89 mm), one-and-half-column (120 mm),
            double-column (183 mm); custom is the width given. A preset sets width. Optional: left out,
            the viewer decides.
        svg_strategy : {"full-vector", "optimized-vector", "hybrid"}
            How an SVG draws its points: full-vector (the default: every point a circle, editable; a
            file grows with the cells drawn), optimized-vector (fewer points, optimized_target_count),
            hybrid (the points as one image, text and axes as vector: small at any cell count).
            Optional: left out, the viewer decides.
        title : str or None
            The figure title while it does not follow the dataset, else null. Optional: left out, the
            viewer decides.
        title_follows_dataset : bool
            The title is the dataset's name (title is then null). Optional: left out, the viewer
            decides.
        width : float
            The figure's width in millimetres; the figure is exactly this wide, at every resolution, and
            as tall as its plots make it (each keeps its view's shape, or its frame's: set a frame to
            change the shape). A width that is not a preset's makes size_preset custom. Optional: left
            out, the viewer decides.
        x_label : str
            The x axis label. Optional: left out, the viewer decides.
        y_label : str
            The y axis label. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('figure-export.set', auto_text_size=auto_text_size, background=background, background_color=background_color, depth_sort=depth_sort, download=download, dpi=dpi, emphasize_selection=emphasize_selection, export_all_views=export_all_views, font_family=font_family, font_sizes=font_sizes, frame=frame, include_axes=include_axes, include_legend=include_legend, legend_position=legend_position, optimized_target_count=optimized_target_count, preview=preview, reference_grid=reference_grid, selection_muted_opacity=selection_muted_opacity, show_orientation=show_orientation, size_preset=size_preset, svg_strategy=svg_strategy, title=title, title_follows_dataset=title_follows_dataset, width=width, x_label=x_label, y_label=y_label)

    def filter_category_all(
        self,
        *,
        field: Mapping[str, Any],
        visible: bool,
    ) -> ActionRequest:
        r"""Show or hide every category.

        Show or hide every category of a categorical column that has open cells, and the cells without a
        category.

        Action ``filter.category.all``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A categorical column of the cells by name (source "obs"). Keys: ``key``, ``source``.
        visible : bool
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('filter.category.all', field=field, visible=visible)

    def filter_category_cell_count(
        self,
        *,
        field: Mapping[str, Any],
        category: str,
        count: int | None,
        seed: int = 0,
    ) -> ActionRequest:
        r"""Keep N cells of a category.

        Show only \`count\` cells of one category, a random sample drawn with \`seed\` (the same seed
        and count show the same cells; a larger count adds cells); null shows all its cells again.

        Action ``filter.category.cell-count``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A categorical column of the cells by name (source "obs"). Keys: ``key``, ``source``.
        category : str
            A category label as the legend shows it.
        count : int or None
            How many of its cells stay shown; null for all of them.
        seed : int
            The sample's seed. Default ``0``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('filter.category.cell-count', field=field, category=category, count=count, seed=seed)

    def filter_category_missing(
        self,
        *,
        field: Mapping[str, Any],
        visible: bool,
    ) -> ActionRequest:
        r"""Show or hide cells without a category.

        Show or hide the cells of a categorical column that have no category (the legend's Missing row).

        Action ``filter.category.missing``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A categorical column of the cells by name (source "obs"). Keys: ``key``, ``source``.
        visible : bool
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('filter.category.missing', field=field, visible=visible)

    def filter_category_visibility(
        self,
        *,
        field: Mapping[str, Any],
        categories: Sequence[str],
        visible: bool,
    ) -> ActionRequest:
        r"""Show or hide categories.

        Show or hide the cells of one or more categories of a categorical column.

        Action ``filter.category.visibility``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A categorical column of the cells by name (source "obs"). Keys: ``key``, ``source``.
        categories : list
            Category labels as the legend shows them.
        visible : bool
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('filter.category.visibility', field=field, categories=categories, visible=visible)

    def filter_continuous_range(
        self,
        *,
        field: Mapping[str, Any],
        min: float,
        max: float,
    ) -> ActionRequest:
        r"""Set a visible range.

        Show only the cells whose value of a continuous field or gene lies in [min, max]; the colors
        follow when they are rescaled to the range.

        Action ``filter.continuous.range``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The continuous field or gene the view in focus is colored by, by name. Keys: ``key``,
            ``source``.
        min : float
            See the action's summary.
        max : float
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('filter.continuous.range', field=field, min=min, max=max)

    def filter_continuous_reset(
        self,
        *,
        field: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Reset a visible range.

        Show every value of a continuous field or gene again.

        Action ``filter.continuous.reset``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The continuous field or gene the view in focus is colored by, by name. Keys: ``key``,
            ``source``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('filter.continuous.reset', field=field)

    def filter_enable(
        self,
        *,
        field: Mapping[str, Any],
        kind: Literal['category', 'continuous', 'outlier'],
        enabled: bool,
    ) -> ActionRequest:
        r"""Turn a filter on or off.

        Turn a field's category, range or outlier filter off (its settings kept) or on again.

        Action ``filter.enable``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The field or gene whose filter it is, by name. Keys: ``key``, ``source``.
        kind : {"category", "continuous", "outlier"}
            See the action's summary.
        enabled : bool
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('filter.enable', field=field, kind=kind, enabled=enabled)

    def filter_outlier(
        self,
        *,
        field: Mapping[str, Any],
        threshold: float,
    ) -> ActionRequest:
        r"""Filter outliers.

        Hide the cells of each category that lie furthest from it in the latent space: those above the
        threshold quantile (1 hides none).

        Action ``filter.outlier``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The categorical column the view in focus is colored by, by name. Keys: ``key``, ``source``.
        threshold : float
            The quantile above which cells are hidden; 1 hides none.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('filter.outlier', field=field, threshold=threshold)

    def filter_remove(
        self,
        *,
        field: Mapping[str, Any],
        kind: Literal['category', 'continuous', 'outlier'],
    ) -> ActionRequest:
        r"""Remove a filter.

        Remove a field's category, range or outlier filter: every cell it hid is shown again.

        Action ``filter.remove``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            The field or gene whose filter it is, by name. Keys: ``key``, ``source``.
        kind : {"category", "continuous", "outlier"}
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('filter.remove', field=field, kind=kind)

    def gene_set_save(
        self,
        *,
        name: str,
        genes: Sequence[str],
        description: str = '',
    ) -> ActionRequest:
        r"""Keep a gene set.

        A named list of genes, kept with the dataset for analyses and sessions; a later set of the same
        name replaces it.

        Action ``gene-set.save``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        name : str
            The set's name.
        genes : list
            Gene names.
        description : str
            What the set is. Default ``''``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('gene-set.save', name=name, genes=genes, description=description)

    def genes_search(
        self,
        *,
        limit: int = 20,
        names: Sequence[str] | Omit = OMIT,
        query: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Search the genes.

        Gene names matching what is typed, best first, as the gene search box ranks them; or, given
        names, which of them the dataset has.

        Action ``genes.search``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        limit : int
            At most this many names (for a query). Default ``20``.
        names : list
            Gene names to look up at once, exact in any case: genes lists those the dataset has, as it
            names them; lacks the rest. Optional: left out, the viewer decides.
        query : str
            What is typed (give it or names). Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('genes.search', limit=limit, names=names, query=query)

    def highlight_category(
        self,
        *,
        field: Mapping[str, Any],
        categories: Sequence[str],
        color: str | Omit = OMIT,
        page: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Highlight categories.

        Add the cells of one or more categories of a field, as the selected view shows them, each
        category as a group of its own, in one step.

        Action ``highlight.category``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A categorical field, by name. Keys: ``key``, ``source``.
        categories : list
            The categories, by their labels as the legend shows them: each becomes a group of the page,
            all in one step. One already on the page is brought up to the cells the view shows now.
        color : str
            The color of the page the group goes on. With \`page\`, the color of the page made (an
            existing page must already have it). Without, the page shown when it has this color, else
            the first page that has it, else a new page of it. Optional: left out, the viewer decides.
        page : str
            The page the group goes on, by name; a page of that name is made when there is none.
            Omitted: the page shown, or with \`color\`, a page of that color. Optional: left out, the
            viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.category', field=field, categories=categories, color=color, page=page)

    def highlight_cells(
        self,
        *,
        cells: CellsLike,
        color: str | Omit = OMIT,
        label: str | Omit = OMIT,
        page: str | Omit = OMIT,
        tool: Literal['annotation', 'lasso', 'proximity', 'knn'] = 'annotation',
    ) -> ActionRequest:
        r"""Highlight cells.

        Add a group of cells to a highlight page: a confirmed selection, or cells an agent or notebook
        names.

        Action ``highlight.cells``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        cells : cells (rows, a boolean mask or obs names)
            Cells by original position, each once, in any order: the 0-based row of obs in the order the
            data was given. A file prepare() wrote holds its rows in priority order, so a row of it read
            as stored is not that position: its obs["cellucid_original_position"] gives the row's (or
            pass obs names).
        color : str
            The color of the page the group goes on. With \`page\`, the color of the page made (an
            existing page must already have it). Without, the page shown when it has this color, else
            the first page that has it, else a new page of it. Optional: left out, the viewer decides.
        label : str
            The group's label; omitted, "Highlight" (the list shows its cell count). A label already on
            the page gets " (2)", " (3)", …. Optional: left out, the viewer decides.
        page : str
            The page the group goes on, by name; a page of that name is made when there is none.
            Omitted: the page shown, or with \`color\`, a page of that color. Optional: left out, the
            viewer decides.
        tool : {"annotation", "lasso", "proximity", "knn"}
            The selection tool that made the cells (annotation for cells named directly). Default
            ``'annotation'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.cells', cells=cells, color=color, label=label, page=page, tool=tool)

    def highlight_clear_all(
        self,
    ) -> ActionRequest:
        r"""Clear every highlight.

        Remove every group of every highlight page; the pages stay.

        Action ``highlight.clear-all``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.clear-all')

    def highlight_group_remove(
        self,
        *,
        group: str,
        page: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Remove a highlight group.

        Remove one group from a highlight page.

        Action ``highlight.group.remove``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        group : str
            A highlight group, by its label as the group list shows it (without the cell count).
        page : str
            The group's page, by name; omitted, the page shown. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.group.remove', group=group, page=page)

    def highlight_group_toggle(
        self,
        *,
        group: str,
        enabled: bool,
        page: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Show or hide a highlight group.

        Turn one group of a page on or off; a group that is off highlights nothing.

        Action ``highlight.group.toggle``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        group : str
            A highlight group, by its label as the group list shows it (without the cell count).
        enabled : bool
            true to show the group, false to hide it.
        page : str
            The group's page, by name; omitted, the page shown. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.group.toggle', group=group, enabled=enabled, page=page)

    def highlight_page_clear(
        self,
        *,
        page: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Clear a highlight page.

        Remove every group of one highlight page; the page stays.

        Action ``highlight.page.clear``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        page : str
            The page, by name; omitted, the page shown. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.page.clear', page=page)

    def highlight_page_combine(
        self,
        *,
        pages: Sequence[str],
        operation: Literal['intersection', 'union'],
        name: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Combine two highlight pages.

        Make a page of the cells in both pages (intersection) or in either (union), and show it.

        Action ``highlight.page.combine``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        pages : list
            The two pages, by name.
        operation : {"intersection", "union"}
            Cells in both pages, or in either.
        name : str
            The new page's name; omitted, "A ∩ B" or "A ∪ B". Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.page.combine', pages=pages, operation=operation, name=name)

    def highlight_page_create(
        self,
        *,
        color: str | Omit = OMIT,
        name: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Make a highlight page.

        Add an empty highlight page and show it.

        Action ``highlight.page.create``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        color : str
            The new page's color; omitted, the next color of the page palette. Optional: left out, the
            viewer decides.
        name : str
            The new page's name; omitted, the first free "Page N". Optional: left out, the viewer
            decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.page.create', color=color, name=name)

    def highlight_page_delete(
        self,
        *,
        pages: Sequence[str],
    ) -> ActionRequest:
        r"""Delete highlight pages.

        Delete one or several highlight pages and their groups, as one step; at least one page stays.

        Action ``highlight.page.delete``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        pages : list
            The pages to delete, by name; at least one page stays (clear the last one instead).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.page.delete', pages=pages)

    def highlight_page_recolor(
        self,
        *,
        page: str,
        color: str,
    ) -> ActionRequest:
        r"""Recolor a highlight page.

        Set the color a highlight page has in its tab, the analyses and the category builder.

        Action ``highlight.page.recolor``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        page : str
            A highlight page, by its name as its tab shows it.
        color : str
            A color as six hex digits, "#e41a1c".

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.page.recolor', page=page, color=color)

    def highlight_page_rename(
        self,
        *,
        page: str,
        name: str,
    ) -> ActionRequest:
        r"""Rename a highlight page.

        Give a highlight page another name, unique among the pages.

        Action ``highlight.page.rename``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        page : str
            A highlight page, by its name as its tab shows it.
        name : str
            The new name.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.page.rename', page=page, name=name)

    def highlight_page_switch(
        self,
        *,
        page: str,
    ) -> ActionRequest:
        r"""Show a highlight page.

        Show another highlight page: its groups are the highlighted cells.

        Action ``highlight.page.switch``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        page : str
            A highlight page, by its name as its tab shows it.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.page.switch', page=page)

    def highlight_range(
        self,
        *,
        field: Mapping[str, Any],
        min: float,
        max: float,
        color: str | Omit = OMIT,
        page: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Highlight a range of values.

        Add the cells whose value of a continuous field or gene lies in a range, as the selected view
        shows them, as a group.

        Action ``highlight.range``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        field : dict
            A continuous field or a gene, by name. Keys: ``key``, ``source``.
        min : float
            The smallest value included.
        max : float
            The largest value included.
        color : str
            The color of the page the group goes on. With \`page\`, the color of the page made (an
            existing page must already have it). Without, the page shown when it has this color, else
            the first page that has it, else a new page of it. Optional: left out, the viewer decides.
        page : str
            The page the group goes on, by name; a page of that name is made when there is none.
            Omitted: the page shown, or with \`color\`, a page of that color. Optional: left out, the
            viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('highlight.range', field=field, min=min, max=max, color=color, page=page)

    def history_anndata_save(
        self,
        *,
        file: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Save AnnData with this state.

        Write a new Cellucid-ready AnnData file (.zarr.zip) of the whole open dataset with the state
        shown applied, as Python and R apply a session: renamed labels, user-defined fields and
        highlight groups as obs columns, gene sets in uns and var, and a record of what was applied;
        never over the dataset's own file.

        Action ``history.anndata.save``. It acts outside the viewer's state and history (a file the
        person receives in the browser, a publication, or what the browser keeps for later visits),
        which Undo does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        file : str
            The new file's name (.zarr.zip is added when it does not end so); "<dataset> with this
            state.zarr.zip" when left out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.anndata.save', file=file)

    def history_bookmark_set(
        self,
        *,
        step: str,
        name: str | None,
    ) -> ActionRequest:
        r"""Bookmark a step.

        Give a step a bookmark name, unique in its tracker (jump to it by name), or remove it with null.
        Bookmark before destructive or heavy work.

        Action ``history.bookmark.set``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        step : str
            A step of the tracker in use, by its id (history_get lists the steps).
        name : str or None
            The bookmark, or null to remove it.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.bookmark.set', step=step, name=name)

    def history_branch_apply(
        self,
        *,
        branch: str,
        name: str | Omit = OMIT,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Apply a branch's steps here.

        Apply a branch's own steps again, by name, on top of the step shown, as a new branch (no values
        are merged; a step that does not apply here is flagged).

        Action ``history.branch.apply``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them).
        name : str
            A name, unique among its kind (compared without case and extra spaces); no length limit.
            Optional: left out, the viewer decides.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.apply', branch=branch, name=name, tracker=tracker)

    def history_branch_archive(
        self,
        *,
        branch: str,
        archived: bool,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Archive a branch.

        Hide a branch from the graph and the pickers, or show it again; it is kept and can be shared.
        Main cannot be archived.

        Action ``history.branch.archive``. It moves the viewer's history, as the person's Undo, Redo or
        a click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them).
        archived : bool
            See the action's summary.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.archive', branch=branch, archived=archived, tracker=tracker)

    def history_branch_cancel(
        self,
    ) -> ActionRequest:
        r"""Cancel a new branch.

        Forget the branch "Branch here" set to start with the next change.

        Action ``history.branch.cancel``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.cancel')

    def history_branch_create(
        self,
        *,
        step: str,
        name: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Branch here.

        Show a step and start a new branch from it with the next change (a branch exists from its first
        step); the name defaults to the next "Branch n".

        Action ``history.branch.create``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        step : str
            A step of the tracker in use, by its id (history_get lists the steps).
        name : str
            A name, unique among its kind (compared without case and extra spaces); no length limit.
            Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.create', step=step, name=name)

    def history_branch_delete(
        self,
        *,
        branch: str,
        tracker: str | Omit = OMIT,
        with_sub_branches: bool = False,
    ) -> ActionRequest:
        r"""Delete a branch.

        Move a branch's own steps (those no other branch starts from) to Recently deleted, where they
        can be restored; with with_sub_branches, the branches that start from it too. Main and locked
        branches cannot be deleted.

        Action ``history.branch.delete``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them).
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.
        with_sub_branches : bool
            Also delete every branch that starts from it. Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.delete', branch=branch, tracker=tracker, with_sub_branches=with_sub_branches)

    def history_branch_lock(
        self,
        *,
        branch: str,
        locked: bool,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Lock a branch.

        Lock a branch (no new step, rename, delete or compaction until unlocked; a change at its tip
        starts a new branch), or unlock it.

        Action ``history.branch.lock``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them).
        locked : bool
            See the action's summary.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.lock', branch=branch, locked=locked, tracker=tracker)

    def history_branch_rename(
        self,
        *,
        branch: str,
        name: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Rename a branch.

        Give a branch another name, unique in its tracker.

        Action ``history.branch.rename``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them).
        name : str
            A name, unique among its kind (compared without case and extra spaces); no length limit.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.rename', branch=branch, name=name, tracker=tracker)

    def history_branch_set_main(
        self,
        *,
        branch: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Make a branch main.

        Make a branch the tracker's main one, what sharing, replay and loading default to.

        Action ``history.branch.set-main``. It moves the viewer's history, as the person's Undo, Redo or
        a click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them).
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.set-main', branch=branch, tracker=tracker)

    def history_branch_switch(
        self,
        *,
        branch: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Switch branch.

        Show a branch's last step and continue on it; an agent's turn continues on it too.

        Action ``history.branch.switch``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them).
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.switch', branch=branch, tracker=tracker)

    def history_branch_without(
        self,
        *,
        step: str,
        name: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Branch without a step.

        Make a new branch from a step's parent that applies the later steps of the path again without
        it, by name; nothing existing changes. The report says what did not apply.

        Action ``history.branch.without``. It moves the viewer's history, as the person's Undo, Redo or
        a click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        step : str
            A step of the tracker in use, by its id (history_get lists the steps).
        name : str
            A name, unique among its kind (compared without case and extra spaces); no length limit.
            Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.branch.without', step=step, name=name)

    def history_clear(
        self,
        *,
        keep_state: bool,
    ) -> ActionRequest:
        r"""Clear the history.

        Start the branch "Fresh start" from the opening step with nothing behind it: keeping the state
        shown (one squashed step from the dataset as opened to it), or from the dataset as opened (the
        branch starts with the next change). Every branch stays, archivable; either way it is kept in
        Recently changed and can be undone there.

        Action ``history.clear``. It moves the viewer's history, as the person's Undo, Redo or a click
        on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        keep_state : bool
            True keeps the state shown as one step; false starts over from the dataset as opened.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.clear', keep_state=keep_state)

    def history_diff(
        self,
        *,
        a: Mapping[str, Any],
        b: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Compare two steps.

        What differs between the state after two steps, part by part (filters, labels, highlights, …);
        the steps may be of two trackers.

        Action ``history.diff``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        a : dict
            The first step. Keys: ``step``, ``tracker`` (optional).
        b : dict
            The second step. Keys: ``step``, ``tracker`` (optional).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.diff', a=a, b=b)

    def history_get(
        self,
        *,
        all_branches: bool = False,
        branch: str | Omit = OMIT,
        cursor: str | Omit = OMIT,
        fields: Sequence[Literal['params', 'notes', 'bookmark']] | Omit = OMIT,
        limit: int = 200,
        since: float | Omit = OMIT,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Read the history.

        The steps taken on this dataset, who took each (person, agent, notebook), branches and the step
        shown.

        Action ``history.get``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        all_branches : bool
            Every branch, not only the path to the step shown. Default ``False``.
        branch : str
            The path to this branch's tip (by its name or id) instead of the path to the step shown; not
            with all_branches. Optional: left out, the viewer decides.
        cursor : str
            Where the page starts: the next of the last page (its last step's id); the first page when
            left out. Optional: left out, the viewer decides.
        fields : list
            What to add to each step: params (each parameter in words), notes, bookmark. Optional: left
            out, the viewer decides.
        limit : int
            The steps of one page (200 when left out: well under what a tool result carries). Default
            ``200``.
        since : float
            Only the steps recorded at or after this time (milliseconds since 1970: a step's "at"), so a
            second read lists what came since. Optional: left out, the viewer decides.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.get', all_branches=all_branches, branch=branch, cursor=cursor, fields=fields, limit=limit, since=since, tracker=tracker)

    def history_group_create(
        self,
        *,
        from_: str,
        to: str,
        name: str,
    ) -> ActionRequest:
        r"""Group steps.

        Name a run of consecutive steps of one branch as a group (a folder; nothing is changed or
        merged). It holds whole turns; a step is in one group at most. Group a multi-step analysis under
        a short name.

        Action ``history.group.create``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        from_ : str
            Its first step.
        to : str
            Its last step.
        name : str
            A name, unique among its kind (compared without case and extra spaces); no length limit.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.group.create', **{'from': from_}, to=to, name=name)

    def history_group_extend(
        self,
        *,
        group: str,
        to: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Add steps to a group.

        Extend a group to a later step of its branch, taking in every step between (whole turns, none of
        another group). A turn's later steps join its group by themselves.

        Action ``history.group.extend``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        group : str
            A group by its name or its id (history_summary lists them).
        to : str
            Its new last step.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.group.extend', group=group, to=to, tracker=tracker)

    def history_group_remove(
        self,
        *,
        group: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Ungroup.

        Remove a group; its steps stay as they are.

        Action ``history.group.remove``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        group : str
            A group by its name or its id (history_summary lists them).
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.group.remove', group=group, tracker=tracker)

    def history_group_rename(
        self,
        *,
        group: str,
        name: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Rename a group.

        Give a group another name, unique in its tracker.

        Action ``history.group.rename``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        group : str
            A group by its name or its id (history_summary lists them).
        name : str
            A name, unique among its kind (compared without case and extra spaces); no length limit.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.group.rename', group=group, name=name, tracker=tracker)

    def history_journal_empty(
        self,
        *,
        entry: str | Omit = OMIT,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Empty Recently deleted.

        Drop what Recently deleted keeps (one entry, or all of them) so its values can be freed; it can
        no longer be restored.

        Action ``history.journal.empty``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        entry : str
            An entry of Recently deleted, by its id (history_summary counts them). Optional: left out,
            the viewer decides.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.journal.empty', entry=entry, tracker=tracker)

    def history_journal_undo(
        self,
        *,
        entry: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Restore from Recently deleted.

        Undo an edit of the history (a delete, a rename, a note, …) from Recently deleted: exactly as it
        was when nothing changed it since, else as new objects ("X (restored)", removed steps applied
        again as a branch). Undoing a rename whose old name was taken since is refused, naming the name
        to use instead.

        Action ``history.journal.undo``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        entry : str
            An entry of Recently deleted, by its id (history_summary counts them).
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.journal.undo', entry=entry, tracker=tracker)

    def history_jump(
        self,
        *,
        at: float | Omit = OMIT,
        bookmark: str | Omit = OMIT,
        step: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Show a step.

        Show the state after any step of the tracker in use, on any branch, as clicking it in the
        timeline does; nothing is undone or lost. Name exactly one of step, bookmark or at.

        Action ``history.jump``. It moves the viewer's history, as the person's Undo, Redo or a click on
        a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        at : float
            A time (milliseconds since 1970): show the last step at or before it on the branch shown.
            Optional: left out, the viewer decides.
        bookmark : str
            A bookmark's name (any case): show the step it marks. Optional: left out, the viewer
            decides.
        step : str
            A step of the tracker in use, by its id (history_get lists the steps). Optional: left out,
            the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.jump', at=at, bookmark=bookmark, step=step)

    def history_methods(
        self,
        *,
        branch: str | Omit = OMIT,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Write the methods paragraph.

        The methods paragraph of a branch (main when left out): the dataset and its citation, the cells
        used, each step that changed the analysis, the tracker's notes as its opening, and the Cellucid
        version.

        Action ``history.methods``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them). Optional: left out, the viewer
            decides.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.methods', branch=branch, tracker=tracker)

    def history_note_add(
        self,
        *,
        target: Mapping[str, Any],
        text: str,
        from_reply: bool = False,
    ) -> ActionRequest:
        r"""Add a note.

        Write a note on a step, a branch, a turn or a tracker. End each turn that changed something with
        a one-line summary note on the turn.

        Action ``history.note.add``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        target : dict or dict or dict or dict or dict
            What the note is on: a step, a branch, a group, a turn, or a tracker itself.
        text : str
            Markdown: emphasis, inline code and http(s) links; shown without HTML.
        from_reply : bool
            Taken from the agent's final reply (the turn summary written for an agent that wrote none).
            Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.note.add', target=target, text=text, from_reply=from_reply)

    def history_note_edit(
        self,
        *,
        note: str,
        text: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Edit a note.

        Change a note's text. An agent edits only its own notes; the person edits any.

        Action ``history.note.edit``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        note : str
            A note, by its id (history_step and history_summary list them).
        text : str
            Markdown: emphasis, inline code and http(s) links; shown without HTML.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.note.edit', note=note, text=text, tracker=tracker)

    def history_note_remove(
        self,
        *,
        note: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Remove a note.

        Remove a note. An agent removes only its own notes; the person removes any.

        Action ``history.note.remove``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        note : str
            A note, by its id (history_step and history_summary list them).
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.note.remove', note=note, tracker=tracker)

    def history_pseudocode(
        self,
        *,
        branch: str | Omit = OMIT,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Write the pseudocode.

        A branch's steps (main when left out) as pseudocode, one line each, as recorded.

        Action ``history.pseudocode``. It reads only and changes nothing; ``request.result()`` is what
        it reads.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them). Optional: left out, the viewer
            decides.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.pseudocode', branch=branch, tracker=tracker)

    def history_recording_pause(
        self,
        *,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Pause recording.

        Pause the history of a tracker: changes still apply but make no step until it is resumed, and an
        agent's change cannot be undone on its own meanwhile.

        Action ``history.recording.pause``. It moves the viewer's history, as the person's Undo, Redo or
        a click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.recording.pause', tracker=tracker)

    def history_recording_resume(
        self,
        *,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Resume recording.

        Record again: what differs from the step shown becomes one squashed step ("While paused: …"), by
        everyone who acted meanwhile; nothing when nothing differs.

        Action ``history.recording.resume``. It moves the viewer's history, as the person's Undo, Redo
        or a click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.recording.resume', tracker=tracker)

    def history_redo(
        self,
    ) -> ActionRequest:
        r"""Redo.

        Show the step after the one shown now on the branch last taken, as the State panel's Redo does.

        Action ``history.redo``. It moves the viewer's history, as the person's Undo, Redo or a click on
        a step does, and records no step; ``request.result()`` is its outcome.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.redo')

    def history_replay(
        self,
        *,
        saved: str,
        branch: str | Omit = OMIT,
        mapping: Mapping[str, Any] | Omit = OMIT,
        name: str | Omit = OMIT,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Replay a saved history.

        Apply a saved state's branch (its main branch when left out) again here, step by step by name,
        in a new tracker linked to it; with a mapping, other names are used for fields, categories and
        genes. The report says what did not match.

        Action ``history.replay``. It moves the viewer's history, as the person's Undo, Redo or a click
        on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        saved : str
            A saved state of this browser (history_saved_list).
        branch : str
            Its branch; its main branch when left out. Optional: left out, the viewer decides.
        mapping : dict
            Names to use in place of the recorded ones, each from its recorded name (the replay report
            offers them: a step's missing names what matched nothing here). A category is named by its
            field's recorded key. Keys: ``categories``, ``fields``, ``genes``. Optional: left out, the
            viewer decides.
        name : str
            A name, unique among its kind (compared without case and extra spaces); no length limit.
            Optional: left out, the viewer decides.
        tracker : str
            Its tracker, by name or id; its tracker shown when left out. Optional: left out, the viewer
            decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.replay', saved=saved, branch=branch, mapping=mapping, name=name, tracker=tracker)

    def history_report(
        self,
        *,
        as_: Literal['file', 'text'] = 'file',
        file: str | Omit = OMIT,
        pictures: Literal['none', 'default'] | Sequence[str] = 'none',
    ) -> ActionRequest:
        r"""Make the printable report.

        The State panel's printable report of the branch shown up to the step shown (the dataset, the
        methods text, every step with who took it, its parameters, what it changed and its notes, then
        the pseudocode), as an .html file the person receives (asked first) or as its text; with
        pictures of the views after the steps chosen (a file only).

        Action ``history.report``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        as_ : {"file", "text"}
            file: an .html file downloaded to the person's computer, which they open, print or save as
            PDF; text: the same HTML document in the result (no pictures). Default ``'file'``.
        file : str
            The file's name; the dataset's name and "history report" when left out. Optional: left out,
            the viewer decides.
        pictures : {"none", "default"} or list
            The steps that get a picture of the views after them (each is shown in turn and drawn, then
            the history returns to the step shown). A file only. Default ``'none'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.report', **{'as': as_}, file=file, pictures=pictures)

    def history_save(
        self,
        *,
        content: Literal['history', 'state', 'steps'] = 'history',
        name: str | Omit = OMIT,
        results: Literal['include', 'recompute'] = 'include',
    ) -> ActionRequest:
        r"""Save the history in this browser.

        Keep the history (every tracker), the current state only, or the steps only as a named saved
        state of this browser (history_saved_list lists them; replay or load it later).

        Action ``history.save``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        content : {"history", "state", "steps"}
            history: the trackers and branches with their values; state: the current state only, without
            its history (one squashed step "Shared state"); steps: the steps only, no cells or values (a
            replay plan for any dataset). Default ``'history'``.
        name : str
            Its name; the dataset's name and the time when left out. Optional: left out, the viewer
            decides.
        results : {"include", "recompute"}
            Analysis results: include them, or leave them out to compute again on open from their
            recorded inputs (smaller; always so for steps). Default ``'include'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.save', content=content, name=name, results=results)

    def history_saved_delete(
        self,
        *,
        saved: str,
    ) -> ActionRequest:
        r"""Delete a saved state.

        Remove one saved state of this browser (an autosave or a named save, by its id from
        history_saved_list) with every value only it held; it cannot be brought back. The history shown
        stays open.

        Action ``history.saved.delete``. It acts outside the viewer's state and history (a file the
        person receives in the browser, a publication, or what the browser keeps for later visits),
        which Undo does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        saved : str
            Its id (history_saved_list).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.saved.delete', saved=saved)

    def history_saved_list(
        self,
    ) -> ActionRequest:
        r"""List saved states.

        The histories kept in this browser (each dataset's autosave and the named saves): name, dataset,
        trackers, steps, size, and whether this version opens them.

        Action ``history.saved.list``. It reads only and changes nothing; ``request.result()`` is what
        it reads.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.saved.list')

    def history_saved_open(
        self,
        *,
        saved: str,
    ) -> ActionRequest:
        r"""Open a saved state.

        Open one saved state of this browser (by its id from history_saved_list) to continue it, as
        Load…'s Open does: its dataset first when another is shown (with the cells it was saved on),
        then its history, which goes on where it left off (an autosave keeps being written; a named save
        stays as it was and the history goes on as a new autosave). The history shown before stays kept
        in this browser while Autosave keeps it (steps it does not keep are asked about first). To apply
        its steps to the dataset shown instead, use history_replay. A whole-file download is asked about
        first.

        Action ``history.saved.open``. It opens another dataset in place of the one the viewer shows, as
        an agent in the viewer does. A notebook's viewer shows the dataset its notebook serves, so the
        viewer refuses it from Python: show another dataset with :func:`cellucid.show_anndata` or
        :func:`cellucid.serve_anndata`.

        Parameters
        ----------
        saved : str
            Its id (history_saved_list).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.saved.open', saved=saved)

    def history_search(
        self,
        *,
        bookmarked: bool | Omit = OMIT,
        by: Literal['person', 'agent', 'python', 'r', 'replay', 'outside'] | Omit = OMIT,
        domain: Literal['cells', 'mapping', 'expression', 'field-overlays', 'gene-sets', 'views', 'active-field', 'filters', 'colors', 'visibility', 'highlights', 'camera', 'camera-path', 'settings', 'analysis', 'layout'] | Omit = OMIT,
        flagged: bool | Omit = OMIT,
        has_note: bool | Omit = OMIT,
        kind: Literal['state', 'view', 'layout'] | Omit = OMIT,
        step_kind: Literal['root', 'action', 'outside', 'squashed'] | Omit = OMIT,
        text: str | Omit = OMIT,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Search the history.

        The steps whose label, parameters, notes or bookmark match, filtered by who took them, their
        kind, the part of the state they changed, notes, bookmarks and flags.

        Action ``history.search``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        bookmarked : bool
            Optional: left out, the viewer decides.
        by : {"person", "agent", "python", "r", "replay", "outside"}
            Who took them. Optional: left out, the viewer decides.
        domain : {"cells", "mapping", "expression", "field-overlays", "gene-sets", "views", "active-field", "filters", "colors", "visibility", "highlights", "camera", "camera-path", "settings", "analysis", "layout"}
            The part of the state they changed. Optional: left out, the viewer decides.
        flagged : bool
            Stale, or no longer applying. Optional: left out, the viewer decides.
        has_note : bool
            Optional: left out, the viewer decides.
        kind : {"state", "view", "layout"}
            What their action changes. Optional: left out, the viewer decides.
        step_kind : {"root", "action", "outside", "squashed"}
            What they are: root, action, outside (changed outside an action) or squashed. Optional: left
            out, the viewer decides.
        text : str
            Words to find (any case). Optional: left out, the viewer decides.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.search', bookmarked=bookmarked, by=by, domain=domain, flagged=flagged, has_note=has_note, kind=kind, step_kind=step_kind, text=text, tracker=tracker)

    def history_share(
        self,
        *,
        chats: bool = False,
        content: Literal['history', 'state', 'steps'] = 'history',
        file: str | Omit = OMIT,
        include: Sequence[Mapping[str, Any]] | Omit = OMIT,
        journal: bool = False,
        notes: bool = True,
        results: Literal['include', 'recompute'] = 'include',
        source: Literal['full', 'name'] = 'full',
        times: Literal['exact', 'dates', 'none'] | Omit = OMIT,
    ) -> ActionRequest:
        r"""Share the history as a file.

        Write a session file to share: the history (trackers and branches chosen; the tracker in use,
        its main and shown branches by default), the current state only, or the steps only; with or
        without notes, conversations, Recently deleted, times and where the dataset opens from.

        Action ``history.share``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        chats : bool
            Conversations kept with the trackers (never with a state). Default ``False``.
        content : {"history", "state", "steps"}
            history: the trackers and branches with their values; state: the current state only, without
            its history (one squashed step "Shared state"); steps: the steps only, no cells or values (a
            replay plan for any dataset). Default ``'history'``.
        file : str
            The file's name; the dataset's name when left out. Optional: left out, the viewer decides.
        include : list
            The trackers to include, each with its branches (its main and shown branches when left out);
            the tracker in use when left out. Optional: left out, the viewer decides.
        journal : bool
            Recently deleted (never with a state). Default ``False``.
        notes : bool
            Notes and bookmarks (never with a state). Default ``True``.
        results : {"include", "recompute"}
            Analysis results: include them, or leave them out to compute again on open from their
            recorded inputs (smaller; always so for steps). Steps only always recompute, whatever is
            given. Default ``'include'``.
        source : {"full", "name"}
            full: the dataset's name and where it opens from; name: its name only. Default ``'full'``.
        times : {"exact", "dates", "none"}
            exact; dates (each rounded down to its day); none. When left out: exact for a history or
            steps, none for a state. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.share', chats=chats, content=content, file=file, include=include, journal=journal, notes=notes, results=results, source=source, times=times)

    def history_squash(
        self,
        *,
        from_: str,
        to: str,
    ) -> ActionRequest:
        r"""Merge steps into one.

        Replace consecutive steps of one branch (from first to last, whole turns; a group's are its
        steps) with one squashed step from the first one's state to the last one's; elsewhere it replays
        as values only. Refused when another branch starts inside; the steps go to Recently deleted. An
        agent squashes only its own turn's steps.

        Action ``history.squash``. It moves the viewer's history, as the person's Undo, Redo or a click
        on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        from_ : str
            The first step to merge.
        to : str
            The last step to merge.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.squash', **{'from': from_}, to=to)

    def history_stale(
        self,
        *,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""List stale results.

        Steps whose result no longer matches its inputs (the cells, the mapping, a highlight changed
        since), each with the step that made it stale and why.

        Action ``history.stale``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.stale', tracker=tracker)

    def history_step(
        self,
        *,
        step: str,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Read a step.

        One step in full: its parameters, pseudocode, who took it and when, what it changed and read,
        its notes and bookmark, and why it is stale.

        Action ``history.step``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        step : str
            A step, by its id (history_get lists the steps).
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.step', step=step, tracker=tracker)

    def history_step_delete(
        self,
        *,
        step: str,
    ) -> ActionRequest:
        r"""Delete a step.

        Move one step to Recently deleted; only the later steps that changed or read what it changed are
        applied again, and one that no longer applies is kept, flagged.

        Action ``history.step.delete``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        step : str
            A step of the tracker in use, by its id (history_get lists the steps).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.step.delete', step=step)

    def history_step_revert(
        self,
        *,
        step: str,
    ) -> ActionRequest:
        r"""Revert a step.

        Add a step that puts back exactly what a step changed, when no later step on the path changed
        the same parts (else use Branch without a step).

        Action ``history.step.revert``. It changes the data's meaning (fields, labels, filters,
        highlights, results), as one step of the viewer's history, which the person can undo;
        ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        step : str
            A step of the tracker in use, by its id (history_get lists the steps).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.step.revert', step=step)

    def history_summary(
        self,
        *,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Summarize the history.

        The history at a glance: trackers, branches (fork step, tip, steps, main, locked, archived),
        bookmarks, turns, stale results and Recently deleted, without the steps; the notes on the
        tracker, each branch, group and turn in full (a shared session's notes are its sender's
        account). Read it before acting on the history.

        Action ``history.summary``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.summary', tracker=tracker)

    def history_tour_start(
        self,
    ) -> ActionRequest:
        r"""Start the guided tour.

        Play the State panel's guided tour of the branch shown, from the opening to the step shown: each
        step is shown in turn with its caption and notes (a camera step is flown to), for a talk or a
        class. It returns once the tour plays; the person pauses, steps or stops it, or
        history_tour_stop does.

        Action ``history.tour.start``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.tour.start')

    def history_tour_stop(
        self,
    ) -> ActionRequest:
        r"""Stop the guided tour.

        Stop the guided tour, as its Stop button does: the history returns to the step shown before it
        began (unless the person made a step meanwhile).

        Action ``history.tour.stop``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.tour.stop')

    def history_tracker_remove(
        self,
        *,
        tracker: str,
    ) -> ActionRequest:
        r"""Remove a tracker.

        Move a tracker to Recently deleted, where it can be restored; when it is the one in use, the one
        used before it is shown first.

        Action ``history.tracker.remove``. It moves the viewer's history, as the person's Undo, Redo or
        a click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        tracker : str
            A tracker by its name or its id (history_summary lists them).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.tracker.remove', tracker=tracker)

    def history_tracker_rename(
        self,
        *,
        tracker: str,
        name: str,
    ) -> ActionRequest:
        r"""Rename a tracker.

        Give a tracker another name, unique among this dataset's trackers.

        Action ``history.tracker.rename``. It moves the viewer's history, as the person's Undo, Redo or
        a click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        tracker : str
            A tracker by its name or its id (history_summary lists them).
        name : str
            A name, unique among its kind (compared without case and extra spaces); no length limit.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.tracker.rename', tracker=tracker, name=name)

    def history_tracker_start(
        self,
        *,
        name: str | Omit = OMIT,
        start: Literal['opened', 'shown'] = 'opened',
    ) -> ActionRequest:
        r"""Start another tracker.

        Start a new tracker, an independent tree of steps on this dataset: from the dataset as opened
        (the default), or from what is shown.

        Action ``history.tracker.start``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        name : str
            A name, unique among its kind (compared without case and extra spaces); no length limit.
            Optional: left out, the viewer decides.
        start : {"opened", "shown"}
            opened: the dataset as it was opened; shown: what is shown now. Default ``'opened'``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.tracker.start', name=name, start=start)

    def history_tracker_use(
        self,
        *,
        tracker: str,
    ) -> ActionRequest:
        r"""Use another tracker.

        Show another tracker's step shown and record into it.

        Action ``history.tracker.use``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        tracker : str
            A tracker by its name or its id (history_summary lists them).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.tracker.use', tracker=tracker)

    def history_turn_undo(
        self,
        *,
        turn: str,
        only_own: bool = False,
    ) -> ActionRequest:
        r"""Undo a turn.

        Take back one agent or notebook turn on the branch shown: show the step before it, or, with
        only_own, delete only the turn's own steps (applying the others again) when someone else's steps
        came in between.

        Action ``history.turn.undo``. It moves the viewer's history, as the person's Undo, Redo or a
        click on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        turn : str
            A turn, by its id: each change an agent makes names its turn in its result (step.turn);
            history_summary lists them all.
        only_own : bool
            Delete only the turn's own steps, keeping the others (to Recently deleted). Default
            ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.turn.undo', turn=turn, only_own=only_own)

    def history_undo(
        self,
    ) -> ActionRequest:
        r"""Undo.

        Show the state before the step shown now, as the State panel's Undo does; the step stays in the
        history, to redo.

        Action ``history.undo``. It moves the viewer's history, as the person's Undo, Redo or a click on
        a step does, and records no step; ``request.result()`` is its outcome.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.undo')

    def history_verify(
        self,
        *,
        branch: str,
        keep: bool = False,
        tracker: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Verify a branch.

        Replay a branch on this dataset into a scratch tracker and compare every part of the state at
        the end with its last step: reproduced exactly, or the first step whose values differ.

        Action ``history.verify``. It moves the viewer's history, as the person's Undo, Redo or a click
        on a step does, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        branch : str
            A branch by its name or its id (history_summary lists them).
        keep : bool
            Keep the scratch tracker afterwards. Default ``False``.
        tracker : str
            A tracker by its name or its id (history_summary lists them); the tracker in use when left
            out. Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('history.verify', branch=branch, keep=keep, tracker=tracker)

    def layout_get(
        self,
    ) -> ActionRequest:
        r"""Read the layout.

        The sidebar's sections (docked or floating, where, open or closed), whether the sidebar shows
        and how wide, and the theme.

        Action ``layout.get``. It reads only and changes nothing; ``request.result()`` is what it reads.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('layout.get')

    def layout_panel_dock(
        self,
        *,
        panel: str,
    ) -> ActionRequest:
        r"""Dock a section.

        Put a floating sidebar section back into the sidebar, where it was.

        Action ``layout.panel.dock``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        panel : str
            A sidebar section by its title (Visualization, Navigation, …) or its id.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('layout.panel.dock', panel=panel)

    def layout_panel_dock_all(
        self,
    ) -> ActionRequest:
        r"""Dock every section.

        Put every floating sidebar section back into the sidebar, each where it was (the sidebar's Dock
        all button).

        Action ``layout.panel.dock-all``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('layout.panel.dock-all')

    def layout_panel_float(
        self,
        *,
        panel: str,
        height: float | Omit = OMIT,
        left: float | Omit = OMIT,
        top: float | Omit = OMIT,
        width: float | Omit = OMIT,
    ) -> ActionRequest:
        r"""Undock a section.

        Tear a sidebar section off into a floating window, at a given place and size or beside the
        sidebar.

        Action ``layout.panel.float``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        panel : str
            A sidebar section by its title (Visualization, Navigation, …) or its id.
        height : float
            Height in CSS pixels (kept at least its header and first row). Optional: left out, the
            viewer decides.
        left : float
            Distance of its left edge from the left of the window, in CSS pixels. Optional: left out,
            the viewer decides.
        top : float
            Distance of its top edge from the top of the window, in CSS pixels. Optional: left out, the
            viewer decides.
        width : float
            Width in CSS pixels (kept at least its smallest working width). Optional: left out, the
            viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('layout.panel.float', panel=panel, height=height, left=left, top=top, width=width)

    def layout_panel_place(
        self,
        *,
        panel: str,
        left: float,
        top: float,
        width: float,
        height: float,
    ) -> ActionRequest:
        r"""Move or resize a floating section.

        Give a floating sidebar section a place and size; it is kept inside the window.

        Action ``layout.panel.place``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        panel : str
            A sidebar section by its title (Visualization, Navigation, …) or its id.
        left : float
            Distance of its left edge from the left of the window, in CSS pixels.
        top : float
            Distance of its top edge from the top of the window, in CSS pixels.
        width : float
            Width in CSS pixels (kept at least its smallest working width).
        height : float
            Height in CSS pixels (kept at least its header and first row).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('layout.panel.place', panel=panel, left=left, top=top, width=width, height=height)

    def layout_section_open(
        self,
        *,
        section: str,
        open: bool,
    ) -> ActionRequest:
        r"""Open or close a section.

        Open or close a section of the sidebar (docked or floating).

        Action ``layout.section.open``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        section : str
            A sidebar section by its title (Visualization, Navigation, …) or its id.
        open : bool
            True opens it, false closes it.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('layout.section.open', section=section, open=open)

    def layout_sidebar_show(
        self,
        *,
        shown: bool,
    ) -> ActionRequest:
        r"""Show or hide the sidebar.

        Show the sidebar beside the plot, or hide it to give the plot the whole window.

        Action ``layout.sidebar.show``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        shown : bool
            True shows the sidebar, false hides it.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('layout.sidebar.show', shown=shown)

    def layout_sidebar_width(
        self,
        *,
        width: float,
    ) -> ActionRequest:
        r"""Resize the sidebar.

        Set the sidebar's width (kept between its narrowest and widest).

        Action ``layout.sidebar.width``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        width : float
            Width in CSS pixels.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('layout.sidebar.width', width=width)

    def layout_theme_set(
        self,
        *,
        theme: Literal['light', 'dark'],
    ) -> ActionRequest:
        r"""Choose the theme.

        Show the app light or dark.

        Action ``layout.theme.set``. It changes the viewer's layout, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        theme : {"light", "dark"}
            See the action's summary.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('layout.theme.set', theme=theme)

    def mapping_set(
        self,
        *,
        embeddings: Mapping[str, Any] | Omit = OMIT,
        expression: str | None | Omit = OMIT,
        gene_names: str | Omit = OMIT,
        highly_variable: str | None | Omit = OMIT,
        knn: str | None | Omit = OMIT,
        latent: str | None | Omit = OMIT,
        log1p: Sequence[Mapping[str, Any]] | Omit = OMIT,
        modality: str | Omit = OMIT,
        vector_fields: Sequence[Mapping[str, Any]] | Omit = OMIT,
    ) -> ActionRequest:
        r"""Change the data mapping.

        Choose which stored arrays are the embeddings, neighbor graph, genes (modality, expression
        matrix, gene names, log1p, highly variable), latent space and vector fields; groups not named
        stay as they are, except the genes groups, which follow a named modality (its defaults) or
        expression matrix (a gene-name or highly variable column its var axis has is kept, else its
        default), as the panel does. Applied as one change: refused, nothing changes.

        Action ``mapping.set``. It changes the data's meaning (fields, labels, filters, highlights,
        results), as one step of the viewer's history, which the person can undo; ``request.result()``
        is its outcome and ``request.step`` the step.

        Parameters
        ----------
        embeddings : dict
            The embeddings views can show, sorted by key then columns, and the default of each
            dimension. Keys: ``defaults``, ``enabled``. Optional: left out, the viewer decides.
        expression : str or None
            "X", "raw", "layer:<name>", or null for no genes. Optional: left out, the viewer decides.
        gene_names : str
            "index", or "column:<var column>". Optional: left out, the viewer decides.
        highly_variable : str or None
            The boolean var column marking highly variable genes, or null. Optional: left out, the
            viewer decides.
        knn : str or None
            The obsp neighbor graph, or null for none. Optional: left out, the viewer decides.
        latent : str or None
            The obsm latent space outlier filters use, or null. Optional: left out, the viewer decides.
        log1p : list
            Each expression matrix that holds log1p values, with its base (null: e). Optional: left out,
            the viewer decides.
        modality : str
            The modality genes are read from (a dataset of several). Optional: left out, the viewer
            decides.
        vector_fields : list
            Vector fields bound to enabled embeddings, sorted by embedding then key. Optional: left out,
            the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('mapping.set', embeddings=embeddings, expression=expression, gene_names=gene_names, highly_variable=highly_variable, knn=knn, latent=latent, log1p=log1p, modality=modality, vector_fields=vector_fields)

    def my_datasets_forget(
        self,
        *,
        name: str,
        saved_states: bool = False,
    ) -> ActionRequest:
        r"""Forget a dataset of My datasets.

        Remove one dataset from My datasets, as its Forget in My datasets › All… does: the entry and,
        for a file or folder of this computer, the browser's handle of it; its saved states only when
        asked. The dataset itself is not touched, and nothing brings the entry back but opening it
        again.

        Action ``my-datasets.forget``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        name : str
            Its name in My datasets (datasets_list).
        saved_states : bool
            Also delete its saved states (histories kept in this browser). Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('my-datasets.forget', name=name, saved_states=saved_states)

    def setting_set(
        self,
        *,
        setting: str,
        value: str | float | bool,
    ) -> ActionRequest:
        r"""Change a setting.

        Set one render, navigation, overlay or tool setting of the sidebar; settings.list names each,
        with its kind, range and choices.

        Action ``setting.set``. It changes how the cells are shown, as one step of the viewer's history,
        which the person can undo; ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        setting : str
            The setting's id, as settings.list gives it.
        value : str or float or bool
            A number for a slider, true or false for a switch, an option value for a choice, #rrggbb for
            a color.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('setting.set', setting=setting, value=value)

    def settings_list(
        self,
        *,
        section: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""List the settings.

        Every setting of the sidebar with its label, section, kind, value, range or choices, and whether
        it can be changed now.

        Action ``settings.list``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        section : str
            Only the settings of the sidebar section with this title (Visualization, Navigation, …).
            Optional: left out, the viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('settings.list', section=section)

    def stored_list(
        self,
    ) -> ActionRequest:
        r"""List what this browser keeps.

        Everything kept in this browser for later visits, as the Stored dialog lists it: per dataset
        (its file's name and host) the converted copies, downloads, downloaded parts and gene indexes
        with their sizes; then the saved states, My datasets, the connection speeds and the decoder
        choices; and how much room the browser gives the site.

        Action ``stored.list``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('stored.list')

    def stored_remove(
        self,
        *,
        connection_speeds: bool = False,
        dataset: str | Omit = OMIT,
        decoder_choices: bool = False,
        everything: bool = False,
        items: Sequence[str] | Omit = OMIT,
        my_datasets: bool = False,
        saved_states: bool = False,
    ) -> ActionRequest:
        r"""Remove what this browser keeps.

        Remove what the Stored dialog removes: kept items by id, all a dataset's kept items, every saved
        state, My datasets (the list and the files' handles), the connection speeds, the decoder
        choices, or everything (Remove all). An item an open tab uses stays until that tab closes.
        Nothing brings it back; the dataset opens again more slowly the first time.

        Action ``stored.remove``. It acts outside the viewer's state and history (a file the person
        receives in the browser, a publication, or what the browser keeps for later visits), which Undo
        does not take back, and records no step; ``request.result()`` is its outcome.

        Parameters
        ----------
        connection_speeds : bool
            The connection speeds remembered per host. Default ``False``.
        dataset : str
            Every kept item of the dataset with this file name (stored_list). Optional: left out, the
            viewer decides.
        decoder_choices : bool
            The decoder choices measured in this browser. Default ``False``.
        everything : bool
            Everything no open tab uses, as Remove all. Default ``False``.
        items : list
            Kept items by id (stored_list). Optional: left out, the viewer decides.
        my_datasets : bool
            My datasets: the list and the handles of its files (my_datasets_forget forgets one). Default
            ``False``.
        saved_states : bool
            Every saved state of this browser (history_saved_delete deletes one). Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('stored.remove', connection_speeds=connection_speeds, dataset=dataset, decoder_choices=decoder_choices, everything=everything, items=items, my_datasets=my_datasets, saved_states=saved_states)

    def vector_field_summary(
        self,
        *,
        by: Mapping[str, Any],
        categories: Sequence[str] | Omit = OMIT,
        neighbors: int = 30,
        vector_field: str | Omit = OMIT,
        view: str | Omit = OMIT,
        visible_only: bool = False,
    ) -> ActionRequest:
        r"""Summarize where cells are heading.

        Where the cells of each category of a field are heading along a vector field (RNA velocity, a
        drift), in the embedding a view shows: per category, how fast (mean speed relative to the mean
        cell), how aligned its arrows are, their mean direction in the view's axes, and the categories
        its arrows lead into (each cell's arrow shared among its nearest cells by how directly it points
        at them, as scVelo's velocity graph does, here in the embedding). Read from the arrows as drawn
        in that embedding, not from the genes: an embedding can bend, merge or reverse a flow.

        Action ``vector-field.summary``. It reads only and changes nothing; ``request.result()`` is what
        it reads.

        Parameters
        ----------
        by : dict
            The categorical field to group the cells by. Keys: ``key``, ``source``.
        categories : list
            Only these categories of \`by\` (else every category with cells); their arrows still lead
            into every category. Optional: left out, the viewer decides.
        neighbors : int
            How many nearest cells (in the embedding) each cell's arrow is shared among. Default ``30``.
        vector_field : str
            The vector field, by name (dataset_summary lists them). Left out: the one the overlay shows,
            else the dataset's default, else the only one bound to the view's embedding; with several
            bound and none of those, name one. Optional: left out, the viewer decides.
        view : str
            The view whose embedding the arrows are read in: "View N"; the selected view when left out.
            Optional: left out, the viewer decides.
        visible_only : bool
            Only the cells the view shows (its filters and hidden categories), as cells and as
            neighbors. Default ``False``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('vector-field.summary', by=by, categories=categories, neighbors=neighbors, vector_field=vector_field, view=view, visible_only=visible_only)

    def view_cameras_link(
        self,
        *,
        linked: bool,
    ) -> ActionRequest:
        r"""Link cameras.

        Let views that show the same embedding share one camera, or give each view its own.

        Action ``view.cameras.link``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        linked : bool
            true to link, false to unlink.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.cameras.link', linked=linked)

    def view_clear(
        self,
    ) -> ActionRequest:
        r"""Remove every kept view.

        Remove every kept view and show the live view alone (as View 1).

        Action ``view.clear``. It changes how the cells are shown, as one step of the viewer's history,
        which the person can undo; ``request.result()`` is its outcome and ``request.step`` the step.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.clear')

    def view_embedding_cycle(
        self,
        *,
        view: str,
        dimension: Literal[1, 2, 3],
    ) -> ActionRequest:
        r"""Next embedding of a dimension.

        The 1, 2 and 3 keys: a view showing another dimension shows this dimension's default; one
        already showing it moves to the next enabled embedding.

        Action ``view.embedding.cycle``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        view : str
            A view: "View N" as the Compare views list numbers its rows (the live view is View 1 while
            it is shown, the kept views follow in the order they were kept).
        dimension : {1, 2, 3}
            The dimension: 1, 2 or 3.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.embedding.cycle', view=view, dimension=dimension)

    def view_embedding_default(
        self,
        *,
        view: str,
        dimension: Literal[1, 2, 3],
    ) -> ActionRequest:
        r"""Show the default embedding.

        Show the dataset's default embedding of one dimension in a view; the view follows that default
        from then on.

        Action ``view.embedding.default``. It changes how the cells are shown, as one step of the
        viewer's history, which the person can undo; ``request.result()`` is its outcome and
        ``request.step`` the step.

        Parameters
        ----------
        view : str
            A view: "View N" as the Compare views list numbers its rows (the live view is View 1 while
            it is shown, the kept views follow in the order they were kept).
        dimension : {1, 2, 3}
            The dimension: 1, 2 or 3.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.embedding.default', view=view, dimension=dimension)

    def view_embedding_set(
        self,
        *,
        view: str,
        embedding: Mapping[str, Any],
    ) -> ActionRequest:
        r"""Show an embedding.

        Show one enabled embedding in a view; the view keeps it when the data mapping changes its
        defaults.

        Action ``view.embedding.set``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        view : str
            A view: "View N" as the Compare views list numbers its rows (the live view is View 1 while
            it is shown, the kept views follow in the order they were kept).
        embedding : dict
            An embedding: the obsm array and the columns of it plotted. Keys: ``columns``, ``key``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.embedding.set', view=view, embedding=embedding)

    def view_focus(
        self,
        *,
        view: str,
    ) -> ActionRequest:
        r"""Select a view.

        Select the view the sidebar edits and the camera controls drive.

        Action ``view.focus``. It changes how the cells are shown, as one step of the viewer's history,
        which the person can undo; ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        view : str
            A view: "View N" as the Compare views list numbers its rows (the live view is View 1 while
            it is shown, the kept views follow in the order they were kept).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.focus', view=view)

    def view_keep(
        self,
        *,
        branch: str | Omit = OMIT,
        name: str | Omit = OMIT,
    ) -> ActionRequest:
        r"""Keep the view.

        Freeze the selected view as a panel of the multiview, with its coloring, filters and camera, and
        select it; with a branch, the view as it is at the last step of that branch of the history (to
        compare branches side by side).

        Action ``view.keep``. It changes how the cells are shown, as one step of the viewer's history,
        which the person can undo; ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        branch : str
            A branch of the tracker in use, by its name: the kept view shows what the selected view
            shows at its last step (coloring, filters, hidden cells, camera). Optional: left out, the
            viewer decides.
        name : str
            The kept view's name (view_rename gives or changes one later). A name for the view, one line
            of text ("Early", "Day 7, treated"): its row and its pane show it beside its number ("2 ·
            Early"), a figure's panel after its letter ("B. Early": the figure letters its panels A, B,
            … in order, so a name needs no letter of its own). Not "View N". Optional: left out, the
            viewer decides.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.keep', branch=branch, name=name)

    def view_layout_set(
        self,
        *,
        layout: Literal['grid', 'single'],
    ) -> ActionRequest:
        r"""Set the view layout.

        Show every view side by side (grid) or the selected view alone (single).

        Action ``view.layout.set``. It changes how the cells are shown, as one step of the viewer's
        history, which the person can undo; ``request.result()`` is its outcome and ``request.step`` the
        step.

        Parameters
        ----------
        layout : {"grid", "single"}
            "grid" or "single".

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.layout.set', layout=layout)

    def view_picture(
        self,
        *,
        height: int = 480,
        legend: bool = True,
        source: Literal['figure', 'screen', 'export'] = 'figure',
        view: str | Omit = OMIT,
        width: int = 640,
    ) -> ActionRequest:
        r"""Take a picture of a view.

        A PNG of a view as the figure export draws it (points, legend), or as the screen shows it now
        (source "screen": the vector field overlay, the connectivity edges and a playing camera path's
        caption included, which a figure cannot draw), or the whole figure exactly as Export would
        download it now with the Figure export settings (source "export": size, frame, background,
        fonts, title, panel titles, legends, every view while Export all views is on: check a figure
        with it before the person downloads it); returned as base64; nothing is downloaded.

        Action ``view.picture``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Parameters
        ----------
        height : int
            The height of the box the plot fits in, in pixels (source "export": of the box the figure is
            scaled into). Default ``480``.
        legend : bool
            Draw the color legend beside the plot (a figure picture; the screen has none). Default
            ``True``.
        source : {"figure", "screen", "export"}
            figure: one view drawn as the figure export draws it (points, legend; no overlay), with
            plain settings (no title, 12 px text, the screen's background). screen: the view's pane as
            the screen shows it now, overlays included (the moving flow caught at this instant), no
            legend; the view must be on screen. export: the figure Export would download now, with every
            Figure export setting (figure_export_set), scaled down into width × height (ask 1200 × 1200
            or more to read a multi-panel figure's titles). Default ``'figure'``.
        view : str
            The view: "View N" as the Compare views list numbers it; the selected view when left out.
            Not with source "export" (it draws what Export would). Optional: left out, the viewer
            decides.
        width : int
            The width of the box the plot fits in, in pixels: the plot keeps the view's shape inside
            width × height, as the figure export frames it (the picture adds its margins and the
            legend). With source "export": the box the whole figure is scaled down into (never up).
            Default ``640``.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.picture', height=height, legend=legend, source=source, view=view, width=width)

    def view_remove(
        self,
        *,
        view: str,
    ) -> ActionRequest:
        r"""Remove a view.

        Remove a view: a kept view is deleted; the live view (View 1) is hidden and every kept view
        stays. The last open view cannot be removed; view_clear shows the live view again.

        Action ``view.remove``. It changes how the cells are shown, as one step of the viewer's history,
        which the person can undo; ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        view : str
            A view: "View N" as the Compare views list numbers its rows (the live view is View 1 while
            it is shown, the kept views follow in the order they were kept).

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.remove', view=view)

    def view_rename(
        self,
        *,
        view: str,
        name: str | None,
    ) -> ActionRequest:
        r"""Name a view.

        Give a view a name of its own, or take it away (null): its row and its pane show it beside its
        number, a figure's panel title after the panel's letter ("B. Early"). The view is still named
        "View N" in actions.

        Action ``view.rename``. It changes how the cells are shown, as one step of the viewer's history,
        which the person can undo; ``request.result()`` is its outcome and ``request.step`` the step.

        Parameters
        ----------
        view : str
            A view: "View N" as the Compare views list numbers its rows (the live view is View 1 while
            it is shown, the kept views follow in the order they were kept).
        name : str or None
            The new name, or null for none.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('view.rename', view=view, name=name)

    def views_state(
        self,
    ) -> ActionRequest:
        r"""Read the views.

        Each open view: its name (what actions name it by), the name it was given (view_rename), its
        title, the representation it shows, what it is colored by, its filters and camera; the layout;
        the highlight pages with their groups; whether the person's system asks for reduced motion.

        Action ``views.state``. It reads only and changes nothing; ``request.result()`` is what it
        reads.

        Returns
        -------
        ActionRequest
            The request, sent at once (or as soon as the viewer connects). In a notebook widget its
            answer arrives once the cell has finished: collect it in a later cell.
        """
        return self.act('views.state')
