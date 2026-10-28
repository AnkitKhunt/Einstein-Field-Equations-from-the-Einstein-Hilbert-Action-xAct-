(* ============================================================
   Einstein field equations from the action principle with xAct
   S = Int d^4x Sqrt[-g] [ R/(2 kap) + L_matter ]
   kap = 8 Pi G (c = 1), signature (-,+,+,+)
   Matter: (i) minimally coupled scalar  (ii) Maxwell field
   Author: Ankit C Khunt
   Requires: xAct (xTensor), Mathematica >= 10
   ============================================================ *)

Quit[]  (* optional: clean kernel; remove if you run inside a session *)

<< xAct`xTensor`
$PrePrint = ScreenDollarIndices;
$CVVerbose = False;

(* ---------- 1. Geometry ---------- *)
DefManifold[M, 4, {a, b, c, d, e, f, i, j}];
DefMetric[-1, g[-a, -b], CD, {";", "\[Del]"}, PrintAs -> "g"];

(* ---------- 2. Fields and constants ---------- *)
DefConstantSymbol[kap, PrintAs -> "\[Kappa]"];
DefScalarFunction[V];
DefTensor[phi[], M, PrintAs -> "\[Phi]"];
DefTensor[A[-a], M];                       (* Maxwell potential *)
DefTensor[Fa[-a, -b], M, Antisymmetric[{-a, -b}], PrintAs -> "F"];
(* F_ab = nabla_a A_b - nabla_b A_a  (we substitute explicitly below) *)
FRule = Fa[-a_, -b_] :> CD[-a][A[-b]] - CD[-b][A[-a]];

(* ---------- 3. Lagrangian densities ---------- *)
sqrtg = Sqrt[-Detg[]];

LGrav   = RicciScalarCD[]/(2 kap);

LScalar = -1/2 g[a, b] CD[-a][phi[]] CD[-b][phi[]] - V[phi[]];

LMaxwell = (-1/4) g[a, c] g[b, d] Fa[-a, -b] Fa[-c, -d] /. FRule;

(* choose the matter sector: LScalar, LMaxwell, or LScalar + LMaxwell *)
LMatter = LScalar;

LTotal = sqrtg (LGrav + LMatter);

(* ---------- 4. Variation w.r.t. inverse metric g^{ab} ---------- *)
(* VarD integrates by parts (drops boundary terms) *)
Simp[expr_] := expr // ContractMetric // ToCanonical // ScreenDollarIndices;

EOMmetric = VarD[g[a, b], CD][LTotal];
EOMmetric = Simp[EOMmetric/sqrtg];

(* delta S = Int sqrt(-g) [ (1/2kap)G_ab - (1/2) T_ab ] delta g^{ab}
   => Einstein equations: *)
EFE = Simp[2 kap EOMmetric];
Print["Field equations (should read G_ab - kap T_ab = 0):"];
EFE

(* ---------- 5. Gravity part alone ---------- *)
EinsteinPart = Simp[VarD[g[a, b], CD][sqrtg LGrav]/sqrtg 2 kap];
Print["Gravity variation:"]; EinsteinPart
(* check against G_ab *)
Simp[EinsteinPart - (RicciCD[-a, -b] - 1/2 RicciScalarCD[] g[-a, -b])]
(* expect 0 *)

(* ---------- 6. Energy-momentum tensor, T_ab = -2/sqrt(-g) dS_m/dg^{ab} ---------- *)
Tmunu = Simp[-2/sqrtg VarD[g[a, b], CD][sqrtg LMatter]];
Print["Energy-momentum tensor:"]; Tmunu

(* ---------- 7. Check: EFE = G_ab - kap T_ab ---------- *)
Simp[EFE - (RicciCD[-a, -b] - 1/2 RicciScalarCD[] g[-a, -b] - kap Tmunu)]
(* expect 0 *)

(* ---------- 8. Matter field equations ---------- *)
EOMphi = Simp[VarD[phi[], CD][LTotal]/sqrtg];
Print["Scalar field equation:"]; EOMphi   (* box phi - V'(phi) = 0 *)

(* Maxwell: vary A_a (use LMaxwell only) *)
EOMA = Simp[VarD[A[a], CD][sqrtg LMaxwell]/sqrtg];
Print["Maxwell equation:"]; EOMA          (* nabla_b F^{ba} = 0 *)

(* ---------- 9. Bianchi identity / conservation check ---------- *)
(* nabla^a (G_ab) = 0 identically *)
Simp[CD[a][RicciCD[-a, -b] - 1/2 RicciScalarCD[] g[-a, -b]]]
(* expect 0 (contracted Bianchi) *)

(* ---------- 10. Easy modifications ---------- *)
(* f(R) gravity: DefScalarFunction[fR]; LGrav = fR[RicciScalarCD[]]/(2 kap);
   Gauss-Bonnet / quadratic terms: LGrav += al RicciScalarCD[]^2 etc.
   Non-minimal coupling: LMatter += xi RicciScalarCD[] phi[]^2 *)
