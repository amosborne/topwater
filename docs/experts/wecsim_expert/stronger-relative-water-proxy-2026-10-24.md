# Stronger case: cached relative-water geometry proxy

Coordinator record of the Capytaine expert's independent read-only cached-data
review. No new field/solver/geometry integration or physical correction.

Inputs: stronger-buoy-20261024 motion3eae21f9..., headingea2794b3..., and accepted
surface cachecbb1480093d3160a8419cdfa7ca5f6b40af4162d9dd2bf6f289f746b03f837b7.
751 display times at0.08s; global surface grid65x65 with0.1875m spacing.
CG elevation is bilinear interpolation; surface gradients use finite differences.
The meanboat/global embedding and native RzRyRx are each applied once.

With water normal n=(-eta_x,-eta_y,1)/norm and n_body=R_body^T*n,
relative heel=atan2(n_body_y,n_body_z) and trim=atan2(-n_body_x,
sqrt(n_body_y^2+n_body_z^2)). Relative vertical offset is body heave-eta(CG).

- Maximum absolute relative heel3.10259deg at51.20s.
- Maximum absolute relative trim2.87693deg at51.68s.
- Maximum absolute relative vertical offset0.0900471m at46.16s.
- Exact full3001sample global-roll peak6.618565deg occurs51.42s. Nearest cached
  time51.44s has globalroll6.607747deg, relativeheel2.69935deg,
  relativetrim-2.51356deg, relative verticaloffset-0.0186329m, etaCG0.622199m.
- Eleven projected full-hull vertex samples depart from the CG tangent plane
  by up to0.233024m (50.56s); departure at51.44s is0.159446m. This is a sampled
  curvature proxy, not a whole-footprint bound.

Existing exact flatwater/freeheave-trim geometry curves give righting-secant /
tangent-stiffness / wetted-set turnover departures:

| Relative flatwater heel | Secant | Tangent | Wet turnover |
| --- | ---: | ---: | ---: |
|1.5deg|3.14%|7.40%|5.38%|
|2.5deg|5.96%|12.83%|8.53%|
|3.0deg|7.3225%|15.3236%|10.0557%|
|3.5deg|8.6386%|17.6805%|11.5506%|
|5.0deg|12.33%|24.03%|15.87%|

These static percentages are context, NOT error bounds for combined dynamic
heel/trim/heave in curved water. The exploratory5%secant/10%tangent labels were
never adopted as hard case gates. The current~3.1deg relative-heel proxy does
not meet those particularly small-departure labels. Wave following is observed,
but negligible geometry changes, wetting, rim/flooding and nonlinear accuracy
have not been established. The defensible status is a provisional fixed-mean
linear reference with moderate known static departures, not proven linearity.

[Existing geometry evidence](../capytaine_expert/geometry-envelope-2026-10-24.md).
No additional study or physical parameter fitting is authorized by this note.
