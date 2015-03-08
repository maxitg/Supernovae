(* ::Package:: *)

(* ::Title:: *)
(*Constraining spacetime variations of nuclear decay rates from light curves of type Ia supernovae*)


(* ::Subtitle:: *)
(*Ivan Karpikov, Maxim Piskunov & Sergey Troitsky*)


(* ::Text:: *)
(*The luminosity of fading type Ia supernovae is governed by radioactive decays of (\[InvisiblePrefixScriptBase]^56)Ni and (\[InvisiblePrefixScriptBase]^56)Co. The decay rates are proportional to the Fermi coupling constant Subscript[G, F] and, therefore, are determined by the vacuum expectation value v of the Brout\[Dash]Englert\[Dash]Higgs field. We use the publicly available SNLS and UNION2.1 sets of light curves of type Ia supernova at various redshifts to constrain possible spacetime variations of the (\[InvisiblePrefixScriptBase]^56)Ni decay rate. The resulting constraint is not very tight; however, it is the only direct bound on the variation of the decay rate for redshifts up to z ~ 1. We discuss potential applications of the result to searches for non-constancy of Subscript[G, F] and v.*)


BeginPackage["Supernovae`"];


Unprotect[ListImport];


ListImport::usage = StringJoin @ {
	"ListImport[\!\(\*StyleBox[\"file\", \"TI\"]\)] imports data from a .list file, returning a complete Wolfram Language version of it."
};


(* ::Section:: *)
(*Implementation*)


Begin["Supernovae`Private`"];


(* ::Subsection:: *)
(*ListImport*)


SyntaxInformation[ListImport] = {"ArgumentsPattern" -> {_}};


ListImport[fileName_] := Module[
	{contents, meta, description, data, dataAssociation},
	contents = Select[Length @ # > 0 &] @ Import[fileName, "Table"];
	meta = Select[StringQ @ #[[1]] && StringMatchQ[#[[1]], "\\@*"] &] @ contents;
	description = StringTake[#, 2 ;; -2] & /@ Most @ Flatten @ Select[StringQ @ #[[1]] && StringMatchQ[#[[1]], "#*"] &] @ contents;
	data = contents[[Position[contents, {"#end"}][[1, 1]] + 1 ;; ]];
	
	dataAssociation = KeyDrop[#, "Filter"] & /@ GroupBy[(Association @ Thread[description -> #] & /@ data), #[["Filter"]] &];

	Association @ Append[StringTake[#[[1]], 2 ;; ] -> Rest @ # & /@ meta, "LightCurves" -> dataAssociation]
]


End[];


Attributes[ListImport] = {ReadProtected};


Protect[ListImport];


EndPackage[];
