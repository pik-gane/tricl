# Study: memory and network context in face-to-face contact dynamics

Scripts for the analyses of the paper in `../paper/` (SocioPatterns high school 2013, Netzschleuder `sp_high_school`).
Build tricl first (`cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build`), then, from the repository root:

```
python3 study/make_data.py                                   # data/ (contacts with gaps <= 40 s)
python3 study/make_data.py --gap 120 --out study/data_gap120 # sensitivity analysis
python3 study/models.py                                      # configs/M0.yaml ... M8.yaml
for m in 0 1 2 4 3 6 5; do python3 study/fit.py $m; done     # results/fit_M*.json (maximum likelihood, held-out days)
python3 study/c1.py                                          # results/fit_C1.json (competing aging model)
STUDY_DATA=data_gap120 python3 study/fit.py 1 --tag _g120    # (likewise for M0, M2, M4, M3, and c1.py --tag _g120)
python3 study/residuals.py                                   # results/residuals.json (time rescaling)
python3 study/closure.py M2                                  # results/closure_M2.json (simple and pair closure vs simulation)
python3 paper/make_figures.py                                # paper/generated/
```

Models (see `models.py`): M0 timetable; M1 + static context; M2 + network context; M4, M6 = M1 + memory on two or four
time scales; M3, M5 = M2 + the same memory; M7, M8 variants with additive triadic closure. C1 = M1 + power-law aging
(non-Markovian, own implementation in `c1.py`). Fits take minutes (M0-M2) to about an hour (M5) on four cores.
