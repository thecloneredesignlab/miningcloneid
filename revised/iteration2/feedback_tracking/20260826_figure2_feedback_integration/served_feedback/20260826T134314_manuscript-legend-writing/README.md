# Served Feedback

Consumer: `manuscript-legend-writing`

Tracking root: `revised/iteration2/feedback_tracking/20260826_figure2_feedback_integration`

Spans served: 1

Existing items served: 0

Read `served_feedback.json`, inspect the raw feedback spans, existing items, and
relevant manuscript artifacts, then fill out `feedback_response_template.json`.

## Rules

- Spans are raw feedback pointers. Read each served span and decide what, if
  anything, it requires in your scope.
- Use `span_signoffs` to close out spans you have inspected. A signed-off span is
  not reopened for the same consumer.
- Use `new_items` only when one or more spans require a consumer-owned work unit
  in your scope. A new item must cite one or more `source_spans`.
- Use `item_signoffs` to record completion, no-change, or not-applicable work on
  an existing item. Put the consumer's completion summary in `response`.
- Use `item_modifications` only to change the durable item record itself. Every
  actual item modification reopens prior item signoffs.
- Put a given item in only one of `item_signoffs` or `item_modifications`.
- Use `deferred_by_run_policy` only when the routed run policy explicitly
  authorizes deferral; cite its blocker fingerprint in `residual_risk`.
- Sign off an item when your scope is complete for the current item state.
- Do not edit `feedback_tracking/*.jsonl` directly.

Submit the completed response with:

```bash
.agents/skills/feedback-manager/scripts/feedback_signoff \
  /Users/4470246/Downloads/miningcloneid/revised/iteration2/feedback_tracking/20260826_figure2_feedback_integration/served_feedback/20260826T134314_manuscript-legend-writing \
  path/to/feedback_response.json
```
