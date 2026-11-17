# Replies to Overleaf review comments (pilot UQ section)

Changes were made in the source markdown (`datacenter_pilot_report.md`) and the
converter, then the LaTeX was regenerated. Campaigns are now named
**aleatoric-1**, **aleatoric-2** and **faultscreen-1** throughout.

| # | Comment (short) | Reply |
|---|---|---|
| 1 | What is dCCT? | Removed the term. The glossary now says the change in CCT caused by the datacenter is always measured against the same operating condition's baseline. |
| 2 | "Recorded, but not the stability test"? | Nadir and ROCOF were recorded but not used in the analysis, so they are removed from the glossary and the text. |
| 3 | Check what the glossary needs; h-scale unused | Glossary trimmed to terms the report uses. Removed h-scale, nadir, ROCOF, scenario, hosting capacity, stage A/B/C and "uncensored pair". |
| 4 | Explain "the cap" in plain language | Rewritten as "Search cap (0.5 s)": the longest fault we test. Protection systems clear real faults well within 0.5 s, so longer faults are not realistic and we do not search beyond it. |
| 5 | "A result (deliverability binds), not an error"? | Rewritten as "Infeasible dispatch": with the committed generators the dispatch cannot serve the load, so the extra load cannot be delivered under that commitment. This is a finding about the system, not a software error. |
| 6 | Describe the main campaigns in words | Glossary entry "Campaigns": two aleatoric campaigns (aleatoric-1, aleatoric-2) vary the operating condition and the datacenter's location and size at one fault; the fault-screening campaign (faultscreen-1) varies the fault location. |
| 7 | What is the "shipped case"? | It was the ACTIVSg200 dynamic case file distributed with GridKit. The text now says so, and names the two TAMU source files we rebuild it from (PowerWorld `.aux` network file and `.dyr` dynamics file). |
| 8 | Add an intro referencing the figure | Added an introductory paragraph to "Building one operating condition" that walks through Figure 4.1. |
| 9 | How are generators ranked? | Units are ranked in merit order (cheapest first by linear cost, nuclear always on), and committing *n* units means switching on the first *n*. Now stated explicitly. |
| 10 | Why a citation here? | Removed. It cited the SUNDIALS solver paper for a statement about voltage support, which did not fit. |
| 11 | Why verbatim `[18, 38]`? | Now plain text: "18 to 38 units". The converter also no longer typesets numbers and ranges as code. |
| 12 | Why verbatim `n_min`, `n_max`? | Rewritten as "a quarter and three quarters of the way from its lower to its upper end". |
| 13 | "Q-limit" wording | Now "reactive power limits (Q limits)". The AC OPF alternative is reduced to one sentence. |
| 14 | Implementation notes should go at the end | All four implementation-note boxes (dispatch, wind, datacenter, Kestrel) are now at the end of their sections. |
| 15 | What is UPS? | Added UPS (uninterruptible power supply) as an acronym, used via `\ac{ups}`. |
| 16 | Sentence hard to parse; drop "folds" | Replaced by: the datacenter mix is P/I/Z = 80/20/0, with one sentence of rationale. The fold and alphaP discussion and the other mixes tried are removed. |
| 17 | Rewrite in plain English (`ok`, `sim_failed`, `LOCALLY_INFEASIBLE`) | Rewritten: "Simulated counts the runs whose dispatch solved, so that a dynamic simulation was attempted. For the remaining runs no feasible dispatch exists: the committed generators cannot serve the load." |
| 18 | Don't use § | All cross-references are now printed as "Section 9.1.x" via `\ref`. |
| 19 | Use acronyms for DCOPF, ACPF | Done with `\ac{dcopf}`, `\ac{acpf}`, `\ac{acopf}`. Only ACOPF existed in `common/acronyms.tex`; DCOPF, ACPF, CCT and UPS have been added. |
| 20 | Display the formula | Greenwood's formula is now a displayed equation. |
| 21 | Use acronym for CCT | Done with `\ac{cct}`, spelled out on first use. |
| 22 | Why a quote environment? | The working hypothesis is now a normal paragraph. |
| 23 | Human-readable campaign names | Renamed to aleatoric-1, aleatoric-2 and faultscreen-1 in the text, tables, captions and figure titles. |
| 24 | Start with human-readable grid names | Headings now read "Texas synthetic grid (ACTIVSg2000)", and the system section "The Illinois synthetic grid (ACTIVSg200)". |
| 25 | Avoid colloquial language | "Where this goes" is now "Next steps". |

## Other simplifications made with this revision

- Sections on GridKit end-user feedback, limitations and implementation are omitted from the report.
- AC OPF alternative reduced to one sentence; alphaP / fold discussion removed.
- Earlier-campaign material and links to working notes removed.
- Parallel-execution storage details condensed.
