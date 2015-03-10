(* ::Package:: *)

(* ::Title:: *)
(*Constraining spacetime variations of nuclear decay rates from light curves of type Ia supernovae*)


(* ::Subtitle:: *)
(*Ivan Karpikov, Maxim Piskunov & Sergey Troitsky*)


(* ::Text:: *)
(*The luminosity of fading type Ia supernovae is governed by radioactive decays of (\[InvisiblePrefixScriptBase]^56)Ni and (\[InvisiblePrefixScriptBase]^56)Co. The decay rates are proportional to the Fermi coupling constant Subscript[G, F] and, therefore, are determined by the vacuum expectation value v of the Brout\[Dash]Englert\[Dash]Higgs field. We use the publicly available SNLS and UNION2.1 sets of light curves of type Ia supernova at various redshifts to constrain possible spacetime variations of the (\[InvisiblePrefixScriptBase]^56)Ni decay rate. The resulting constraint is not very tight; however, it is the only direct bound on the variation of the decay rate for redshifts up to z ~ 1. We discuss potential applications of the result to searches for non-constancy of Subscript[G, F] and v.*)


BeginPackage["Supernovae`"];
Needs["ErrorBarPlots`"];


Unprotect[ListImport, MaxFluxDate, ExponentialDecaySubset, NormalizeLightCurve, NormalizedExponentialDecayFit, NormalizedExponentialDecayPlot, SupernovaPlot];


ListImport::usage = StringJoin @ {
	"ListImport[\!\(\*StyleBox[\"file\", \"TI\"]\)] imports data from a .list file, returning a complete Wolfram Language version of it."
};


MaxFluxDate::usage = StringJoin @ {
	"MaxFluxDate[\!\(\*StyleBox[\"lightCurve\", \"TI\"]\)] yields the date of maximal flux in lightCurve."
};


ExponentialDecaySubset::usage = StringJoin @ {
	"ExponentialDecayInterval[\!\(\*StyleBox[\"lightCurve\", \"TI\"]\)] yields a subset of lightCurve during which flux exponentially decays after the global maximum."
};


NormalizeLightCurve::usage = StringJoin @ {
	"NormalizeLightCurve[\!\(\*StyleBox[\"lightCurve\", \"TI\"]\)] yields a modification of lightCurve, in which fluxes are replaced with log10 of fluxes, ",
	"and dates are normalized, so that the date of maximal flux is zero. \n",
	"NormalizeLightCurve[\!\(\*StyleBox[\"lightCurve\", \"TI\"], StyleBox[\"zeroDate\", \"TI\"]\)] normalizes the dates, so that zeroDate becomes zero."
};


NormalizedExponentialDecayFit::usage = StringJoin @ {
	"NormalizedExponentialDecayFit[\!\(\*StyleBox[\"lightCurve\", \"TI\"]\)] yields a linear fit model of NormalizeLightCurve[\!\(\*StyleBox[\"lightCurve\", \"TI\"]\)]."
};


NormalizedExponentialDecayPlot::usage = StringJoin @ {
	"NormalizedExponentialDecayPlot[\!\(\*StyleBox[\"lightCurve\", \"TI\"]\)] yields a log plot of lightCurve including fit function."
};


SupernovaPlot::usage = StringJoin @ {
	"SupernovaPlot[\!\(\*StyleBox[\"Supernovae\", \"TI\"]\)] yields a log plot of Supernovae light curves for all filters including fit functions."
};


(* ::Section:: *)
(*Implementation*)


Begin["Supernovae`Private`"];


(* ::Subsection:: *)
(*ListImport*)


SyntaxInformation[ListImport] = {"ArgumentsPattern" -> {_}};


ListImport[fileName_String] /; FileExistsQ[fileName] := Module[
	{contents, meta, description, data, dataAssociation},
	contents = Select[Length @ # > 0 &] @ Import[fileName, "Table"];
	meta = Select[StringQ @ #[[1]] && StringMatchQ[#[[1]], "\\@*"] &] @ contents;
	description = StringTake[#, 2 ;; -2] & /@ Most @ Flatten @ Select[StringQ @ #[[1]] && StringMatchQ[#[[1]], "#*"] &] @ contents;
	data = contents[[Position[contents, {"#end"}][[1, 1]] + 1 ;; ]];
	
	dataAssociation = KeyDrop[#, "Filter"] & /@ GroupBy[(Association @ Thread[description -> #] & /@ data), #[["Filter"]] &];

	Association @ Append[StringTake[#[[1]], 2 ;; ] -> Rest @ # & /@ meta, "LightCurves" -> dataAssociation]
]


(* ::Subsection:: *)
(*$LightCurveSortedByTime*)


(* ::Text:: *)
(*Discussion! Negative fluxes are dropped.*)


$LightCurveSortedByTime[lightCurve_List] := Select[#["Flux"] > 0 &] @ SortBy[lightCurve, #["Date"]&];


(* ::Subsection:: *)
(*MaxFluxDate*)


SyntaxInformation[MaxFluxDate] = {"ArgumentsPattern" -> {_}};


MaxFluxDate[lightCurve_List] := Module[
	{lightCurveSortedByTime, maxFluxDate, lightCurveAfterPeak, logFlux, weights},
	lightCurveSortedByTime = $LightCurveSortedByTime @ lightCurve;
	lightCurveSortedByTime[[Position[lightCurveSortedByTime, Max[lightCurveSortedByTime[[All, "Flux"]]]][[-1, 1]], "Date"]]
]


(* ::Subsection:: *)
(*ExponentialDecayRange*)


SyntaxInformation[ExponentialDecaySubset] = {"ArgumentsPattern" -> {{_}, OptionsPattern[]}};


Options[ExponentialDecaySubset] = {"MaxDurationAfterPeak" -> 30, "MinAdjustedRSquared" -> .95, "FitSelectionMethod" -> "NormalizedChiSquared"};


ExponentialDecaySubset[lightCurve_List, OptionsPattern[]] := Module[
	{maxFluxDate, lightCurveAfterPeak, logFlux, weights},
	maxFluxDate = MaxFluxDate @ lightCurve;
	lightCurveAfterPeak = Select[maxFluxDate < #["Date"] < maxFluxDate + OptionValue["MaxDurationAfterPeak"] &] @ $LightCurveSortedByTime @ lightCurve;
	logFlux = {#[[1]], Log10 @ #[[2]]} & /@ lightCurveAfterPeak;
	weights = #[[3]]/(Log[10] #[[2]]) & /@ lightCurveAfterPeak;
	If[Head @ # === Missing, {}, lightCurveAfterPeak[[#[[1]] ;; #[[2]]]]] & @
		If[OptionValue["FitSelectionMethod"] == "AdjustedRSquared",
			SelectFirst[LinearModelFit[logFlux[[#[[1]] ;; #[[2]]]], t, t, Weights -> weights[[#[[1]] ;; #[[2]]]]]["AdjustedRSquared"] >= OptionValue["MinAdjustedRSquared"] &] @*
			SortBy[#[[1]] - #[[2]] &],
			If[# == {}, Missing[], First @ #] & @*
			SortBy[Total[LinearModelFit[logFlux[[#[[1]] ;; #[[2]]]], t, t, Weights -> weights[[#[[1]] ;; #[[2]]]]]["StandardizedResiduals"]^2]/(#[[2]] - #[[1]] - 1) &]
		] @
		Select[#[[2]] - #[[1]] > 1 &] @
		Subsets[Range @ Length @ logFlux, {2}]
]


(* ::Subsection:: *)
(*NormalizeLightCurve*)


SyntaxInformation[NormalizeLightCurve] = {"ArgumentsPattern" -> {__}};


NormalizeLightCurve[lightCurve_, zeroDate_] := {#[["Date"]] - zeroDate, Log10 @ #[["Flux"]], #[["Fluxerr"]]/(Log[10] #[["Flux"]])} & /@ lightCurve


NormalizeLightCurve[lightCurve_] := NormalizeLightCurve[lightCurve, MaxFluxDate @ lightCurve]


(* ::Subsection:: *)
(*NormalizedExponentialDecayFit*)


SyntaxInformation[NormalizedExponentialDecayFit] = {"ArgumentsPattern" -> {{_}, OptionsPattern[]}};


Options[NormalizedExponentialDecayFit] = {"MaxDurationAfterPeak" -> 30, "MinAdjustedRSquared" -> .95};


NormalizedExponentialDecayFit[lightCurve_List, options___] := Module[
	{normalizedLightCurve},
	normalizedLightCurve = NormalizeLightCurve[ExponentialDecaySubset[$LightCurveSortedByTime @ #, options], MaxFluxDate @ #] & @ lightCurve;
	If[Length @ normalizedLightCurve == 0, Missing["No good fit."], LinearModelFit[normalizedLightCurve[[All, {1, 2}]], t, t, Weights -> normalizedLightCurve[[All, 3]]]]
]


(* ::Subsection:: *)
(*NormalizedExponentialDecayPlot*)


SyntaxInformation[NormalizedExponentialDecayPlot] = {"ArgumentsPattern" -> {{_}, OptionsPattern[]}};


Options[NormalizedExponentialDecayPlot] = {
	"MaxDurationAfterPeak" -> 30,
	"MinAdjustedRSquared" -> .95,
	"FitSelectionMethod" -> "NormalizedChiSquared",
	"ExcludedPlotStyle" -> Gray,
	"FitPlotStyle" -> ColorData[97, 4],
	"PlotMarkers" -> {Automatic, Medium},
	"ImageSize" -> 300,
	"FrameLabel" -> {"Time, days", "Log10[flux, \!\(\*SuperscriptBox[\(sec\), \(-1\)]\)]"},
	"PlotRange" -> All,
	"PlotLabel" -> ""
};


NormalizedExponentialDecayPlot[lightCurve_List, OptionsPattern[]] := Module[
	{maxFluxDate, normalizedLightCurve, normalizedExponentialDecay, fitFunction},
	maxFluxDate = MaxFluxDate @ lightCurve;
	normalizedLightCurve = NormalizeLightCurve[$LightCurveSortedByTime @ lightCurve, maxFluxDate];
	normalizedExponentialDecay = NormalizeLightCurve[
		ExponentialDecaySubset[
			lightCurve,
			"MaxDurationAfterPeak" -> OptionValue["MaxDurationAfterPeak"],
			"MinAdjustedRSquared" -> OptionValue["MinAdjustedRSquared"],
			"FitSelectionMethod" -> OptionValue["FitSelectionMethod"]
		],
		maxFluxDate
	];
	fitFunction = NormalizedExponentialDecayFit[
		lightCurve,
		"MaxDurationAfterPeak" -> OptionValue["MaxDurationAfterPeak"],
		"MinAdjustedRSquared" -> OptionValue["MinAdjustedRSquared"],
		"FitSelectionMethod" -> OptionValue["FitSelectionMethod"]
	];
	
	Show[
		ErrorListPlot[{
			{{#[[1]], #[[2]]}, ErrorBar @ #[[3]]} & /@ Complement[normalizedLightCurve, normalizedExponentialDecay],
			{{#[[1]], #[[2]]}, ErrorBar @ #[[3]]} & /@ normalizedExponentialDecay /. {} -> {Missing[]}
			},
			PlotStyle -> {OptionValue["ExcludedPlotStyle"], OptionValue["FitPlotStyle"]},
			PlotMarkers -> OptionValue["PlotMarkers"],
			Method -> {"OptimizePlotMarkers" -> False}
		],
		If[Length @ normalizedExponentialDecay > 0, Plot[fitFunction[t], {t, normalizedExponentialDecay[[1, 1]], normalizedExponentialDecay[[-1, 1]]}, PlotStyle -> OptionValue["FitPlotStyle"]], {}],
		ImageSize -> OptionValue["ImageSize"],
		Axes -> False,
		Frame -> True,
		FrameLabel -> OptionValue["FrameLabel"],
		PlotRange -> OptionValue["PlotRange"],
		PlotLabel -> OptionValue["PlotLabel"]
	]
]


(* ::Subsection:: *)
(*SupernovaPlot*)


End[];


Attributes[ListImport] = {ReadProtected};
Attributes[MaxFluxDate] = {ReadProtected};
Attributes[ExponentialDecaySubset] = {ReadProtected};
Attributes[NormalizeLightCurve] = {ReadProtected};
Attributes[NormalizedExponentialDecayFit] = {ReadProtected};
Attributes[NormalizedExponentialDecayPlot] = {ReadProtected};
Attributes[SupernovaePlot] = {ReadProtected};


Protect[ListImport, MaxFluxDate, ExponentialDecaySubset, NormalizeLightCurve, NormalizedExponentialDecayFit, NormalizedExponentialDecayPlot, SupernovaPlot];


EndPackage[];
