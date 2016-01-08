(* ::Package:: *)

(* :Title: shasta *)

(* :Context: shasta` *)

(* :Author: Greg Jackson (greg@wam.co.za) *)

(* :Summary:

  SHASTA = a SHarp And Smooth Transport Algorithm. It is designed to handle
  shock waves in transport equations. Here the package is intended for use
  in Heavy Ion Collisions. In the hydrodynamic limit, L >> mfp, the energy-
  momentum tensor is described by equilibrium parameters e & p. Conservation
  laws \!\(
\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]
\*SuperscriptBox[\(T\), \(\[Mu]\[Nu]\)]\) = 0 describe transport of the components.

  Historically, relativistic transport equations where numerically studied in
  the 1970s (with application to nuclear reactions). The first paper describing
  the SHASTA algorithm is [Boris & Book, J. Comp. Phys. 11, 38-69 (1973)]. 
  Shortly afterwards, [Zalesak, J. Comp. Phys. 31, 335-362 (1979)] generalised 
  the scheme to multidimensional transport equations.

*)

(* :Package Version: 1.0 *)

(* :Mathematica Version: 10.0 *)

(* :History:
        
  V1.0 2015/04/12:
    Includes a 1-d & 2-d flux limiter
          
*)

(* :Keywords:
  flux-limit, transport, shock
*)

(* :Limitations:  
  
*)

BeginPackage["shasta`"]

Bare1D::usage=
  "Bare1D[U,V,\[Lambda]] updates by one timestep \[CapitalDelta]t the solution to the (1d) transport eqaution."

Bare2D::usage=
  "Bare2D[U,V,W,\[Lambda]x,\[Lambda]y] updates by one timestep \[CapitalDelta]t the solution to the (2-d) transport eqaution."

FluxL1D::usage=
  "FluxL1D[U] applies the SHASTA correction to a (1-d) list U. Intended for use with Bare1D[...]."

FluxL2D::usage=
  "FluxL2D[U] applies the SHASTA correction to a (2-d) array U. Intended for use with Bare2D[...]."

AntiD::usage=
  "AntiD[U] returns the antidiffusive fluxes for (1-d) list U."

MaxMin::usage=
  "MaxMin[x,y,z] = Sign[y]Max[0,Min[Sign[y]x,(1/8)Abs[y],Sign[y]z]], where x,y,z are lists."


Begin[ "Private`"]

protected=Unprotect[Bare1D, Bare2D, AntiD, Q, FluxL1D, FluxL2D]

(*Clear[protected];*)


(* Here I define many helper functions *)

Extend = Join[{#[[1]]},#,{#[[-1]]}]&;

Grow = (#//Extend//Transpose//Extend//Transpose)&;

Shrink = #[[2;;-2,2;;-2]]&; (* remove ghosts *)

Seq = 1+UnitStep[#1];;-1-UnitStep[-#1]&;

DiffE = Differences[#//Extend]&;

Mid = #[[2;;-2]]&;

MaxMin[x_,y_,z_] := Sign[y]Max[0,Min[Sign[y]x,(1/8)Abs[y],Sign[y]z]];

\[Theta]Min[x_,y_,z_] = UnitStep[x]Min[1,y,z];

Lmax[x_,y_,z_] := Max[x,y,z];
Lmax[x_,y_] := Max[x,y];
Lmin[x_,y_,z_] := Min[x,y,z];
Lmin[x_,y_] := Min[x,y];

SetAttributes[#,Listable]&/@{MaxMin,\[Theta]Min,Lmax,Lmin};

Lf[l_,u_] := (l[ #[[3;;]],#[[2;;-2]],#[[;;-3]] ]&/@u)//Mid;

MmM[l_,\[Sigma]_] := (Max[0,#]&/@l[[ Seq[-\[Sigma]] ]])-( Min[0,#]&/@l[[ Seq[\[Sigma]] ]]);
in = MmM[#[[2;;-2]],+1]&;
out = MmM[#[[2;;-2]],-1]&;


Q[v_,\[Lambda]_,\[Sigma]_] := ((1/2 \[Lambda]^-1-\[Sigma] v[[Seq[-\[Sigma]]]])/(\[Lambda]^-1+Differences[v]))[[Seq[\[Sigma]]]];

Bare1D[u_,v_,\[Lambda]_] := Module[
  {
    \[Delta] = Differences[u//Extend], vE=v//Extend,
    $Qsub = Subscript[(q_), \[Sigma]_:IntegerQ][arg__]:>q[arg,\[Sigma]]
  },
  ReleaseHold[(
    (Hold[# Subscript[Q, #][vE,\[Lambda]]^2 \[Delta][[Seq[#]]]/2 + Subscript[Q, #][vE,\[Lambda]]u]&/@{+1,-1})/.$Qsub
  )//Total]
];

Bare2D[u_,v_,w_,\[Lambda]x_,\[Lambda]y_] := 
  ((Bare1D[
    (u//Transpose)[[#]], (v//Transpose)[[#]], \[Lambda]x
  ]&/@Range[1,Length[u//Transpose]])//Transpose)+
  ((Bare1D[
    (u[[#]]), (w[[#]]), \[Lambda]y
  ]&/@Range[1,Length[u]]))-u;


AntiD[u_] := MaxMin[{0}~Join~#[[;;-2]],#,#[[2;;]]~Join~{0}]&[DiffE[u]];

FluxL1D[u_] := u - Differences[AntiD[u]];

FluxL2D[Ut_] := Module[
  {
    tuu, \[ScriptCapitalA], Aio, F, A, \[ScriptCapitalU]mm, fnc, $Ssub
  },
  $Ssub=fnc[\[Sigma]_:IntegerQ]:>#1 \[Theta]Min@@{ \[Sigma] #1,#2[[1, Seq[+\[Sigma]] ]],#2[[2, Seq[-\[Sigma]] ]] }&;
  tuu = {Transpose[#],#}&[Ut//Grow]; \[ScriptCapitalA] = (AntiD/@#)&/@ tuu;
  \[ScriptCapitalU]mm = Apply@@{
    #, Transpose[{Transpose[#2],#1}&@@(Lf[fnc,#]&/@tuu/.fnc->#),{3,2,1}],{2}
  }&/@{Lmax, Lmin};
  Aio = Mid[#2]+(Mid[#1]//Transpose)&@@((fnc/@#)&/@\[ScriptCapitalA]/.fnc->#)&/@{in,out};
  F = Grow/@((\[ScriptCapitalU]mm-{Ut,Ut}){+1,-1}/(Aio/.x_/;x==0->1));
  A = ((fnc[#]&/@{+1,-1}//Total)/.$Ssub@@ #&/@Transpose[
    { Mid[#//Transpose]&/@\[ScriptCapitalA],{F,Transpose/@F} }
  ]);
  (*(Ut) - (#2) - Transpose[#1] & @@ (Shrink/@(DiffE/@\[ScriptCapitalA][[1]]))*)
  (* Dimensions[Total[Shrink/@{ DiffE/@Transpose[A[[1]]],DiffE/@Transpose[A[[2]]] }]]*)
  (Ut) - Total[Shrink/@{ DiffE/@Transpose[A[[2]]],Transpose[DiffE/@Transpose[A[[1]]]] }]
];


Protect[ Evaluate[protected] ];

End[];
 
EndPackage[];
