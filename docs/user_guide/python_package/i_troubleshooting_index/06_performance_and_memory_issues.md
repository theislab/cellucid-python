# Performance and memory issues

First identify the phase that is slow:

1. Python preparation and validation
2. writing or reading export artifacts
3. network or disk transfer
4. browser parsing and GPU upload
5. interaction, filtering, analysis, or figure export

Record dataset shape, sparse or dense representation, artifact sizes, browser,
viewport, number of visible views, render mode, and the exact operation being
timed.

## Python memory grows during preparation

Preserve sparse matrices when the API accepts them, avoid unnecessary dense
copies, and verify expression orientation before export.

## The browser becomes slow

Use the in-app Performance Benchmark and browser performance tools to distinguish
frame-rate pressure from long CPU tasks or slow requests. Additional views and
large pixel dimensions add independent rendering work.

## Figure export is slow or large

Record format, dimensions, DPI, visible point count, and whether all views are
included. Test the smallest scientifically sufficient output first.

See {doc}`../b_concepts_mental_models/07_performance_mental_model_and_scaling`,
{doc}`../c_data_preparation_api/10_performance_tuning_guide_prepare_export`, and
{doc}`../d_viewing_apis/14_performance_scaling_and_lazy_loading`.
