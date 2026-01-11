// MMGBSA.cpp - MM/GBSA Implementation
#include "MMGBSA.h"
#include "PDBViewer.h"
#include "Components/StaticMeshComponent.h"
#include "Kismet/GameplayStatics.h"
#include "Misc/FileHelper.h"
#include "Misc/Paths.h"

AMMGBSA::AMMGBSA()
{
    PrimaryActorTick.bCanEverTick = false;
}

void AMMGBSA::BeginPlay()
{
    Super::BeginPlay();
    InitializeElementDatabase();

    // Auto-find viewer
    TArray<AActor *> Found;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), APDBViewer::StaticClass(), Found);
    if (Found.Num() > 0)
    {
        ViewerReference = Cast<APDBViewer>(Found[0]);
        // If a viewer exists, load current atom data so calculations can proceed immediately
        LoadAtomsFromViewer();
        // Bind to ligand loaded event so we refresh cache when the viewer updates
        ViewerReference->OnLigandsLoaded.AddDynamic(this, &AMMGBSA::OnViewerLigandsLoaded);
    }
}

void AMMGBSA::InitializeElementDatabase()
{
    ElementDatabase.Empty();

    auto Add = [this](const FString &Sym, float Rad, float Eps, float Chg, float GBScale)
    {
        FElementData Data;
        Data.VDWRadius = Rad;
        Data.VDWEpsilon = Eps;
        Data.TypicalCharge = Chg;
        Data.GBRadiusScale = GBScale;
        ElementDatabase.Add(Sym, Data);
    };

    // Radii in Ångstroms (not scaled), epsilon in kcal/mol
    Add(TEXT("H"), 1.20f, 0.0157f, 0.0f, 0.85f);
    Add(TEXT("C"), 1.70f, 0.0860f, 0.0f, 1.00f);
    Add(TEXT("N"), 1.55f, 0.1700f, -0.3f, 1.05f);
    Add(TEXT("O"), 1.52f, 0.2100f, -0.4f, 1.10f);
    Add(TEXT("S"), 1.80f, 0.2500f, -0.2f, 1.15f);
    Add(TEXT("P"), 1.80f, 0.2000f, 0.5f, 1.10f);
    Add(TEXT("F"), 1.47f, 0.0610f, -0.2f, 1.00f);
    Add(TEXT("CL"), 1.75f, 0.2650f, -0.1f, 1.05f);
    Add(TEXT("BR"), 1.85f, 0.3200f, -0.1f, 1.10f);
    Add(TEXT("I"), 1.98f, 0.4000f, -0.1f, 1.15f);

    // Metal ions (common in protein structures)
    Add(TEXT("ZN"), 1.39f, 0.2500f, 2.0f, 1.20f); // Zinc(II)
    Add(TEXT("Z"), 1.39f, 0.2500f, 2.0f, 1.20f);  // Zinc alias
    Add(TEXT("MG"), 1.73f, 0.8750f, 2.0f, 1.18f); // Magnesium(II)
    Add(TEXT("M"), 1.73f, 0.8750f, 2.0f, 1.18f);  // Magnesium alias
    Add(TEXT("CA"), 2.31f, 0.4500f, 2.0f, 1.22f); // Calcium(II)
    Add(TEXT("FE"), 1.47f, 0.0500f, 2.0f, 1.20f); // Iron(II/III)
    Add(TEXT("MN"), 1.61f, 0.0500f, 2.0f, 1.20f); // Manganese(II)
    Add(TEXT("CU"), 1.40f, 0.0500f, 2.0f, 1.20f); // Copper(II)
    Add(TEXT("NA"), 2.27f, 0.0028f, 1.0f, 1.15f); // Sodium(I)
    Add(TEXT("K"), 2.75f, 0.0003f, 1.0f, 1.25f);  // Potassium(I)
}

void AMMGBSA::InitializeFromViewer(APDBViewer *Viewer)
{
    if (!Viewer)
    {
        UE_LOG(LogTemp, Error, TEXT("MMGBSA: Null viewer"));
        return;
    }

    ViewerReference = Viewer;
    LoadAtomsFromViewer();

    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Loaded %d receptor atoms, %d ligands"),
           ReceptorAtoms.Num(), LigandAtoms.Num());
}

void AMMGBSA::LoadAtomsFromViewer()
{
    if (!ViewerReference)
        return;

    ReceptorAtoms.Empty();
    LigandAtoms.Empty();
    CachedResults.Empty();

    if (ViewerReference)
    {
        ViewerReference->ClearOverlapMarkers();
    }
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Loading atoms from viewer - ligands: %d"), ViewerReference->LigandMap.Num());

    // Load receptor atoms (from ResidueMap - these are ATOM entries)
    {
        TSet<FString> UnknownElementsReceptor;
        for (const auto &Pair : ViewerReference->ResidueMap)
        {
            const FResidueInfo *ResInfo = Pair.Value;
            if (!ResInfo || !ResInfo->bIsVisible)
                continue;

            for (int32 i = 0; i < ResInfo->AtomPositions.Num(); ++i)
            {
                FMMAtom Atom;
                Atom.Position = ResInfo->AtomPositions[i];
                Atom.bIsReceptor = true;

                if (ResInfo->AtomElements.IsValidIndex(i) && !ResInfo->AtomElements[i].IsEmpty())
                    Atom.Element = ResInfo->AtomElements[i].ToUpper();
                else
                    Atom.Element = TEXT("C");

                const FElementData *Data = ElementDatabase.Find(Atom.Element);
                if (Data)
                {
                    Atom.Radius = Data->VDWRadius * 100.0f;
                    Atom.Charge = Data->TypicalCharge * ChargeScaling; // Apply charge scaling to receptor atoms
                }
                else
                {
                    if (!UnknownElementsReceptor.Contains(Atom.Element))
                    {
                        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Unknown element '%s' encountered when loading receptor residue %s. Defaulting to C."), *Atom.Element, *Pair.Key);
                        UnknownElementsReceptor.Add(Atom.Element);
                    }
                    Atom.Element = TEXT("C");
                    const FElementData *Def = ElementDatabase.Find(Atom.Element);
                    if (Def)
                    {
                        Atom.Radius = Def->VDWRadius * 100.0f;
                        Atom.Charge = Def->TypicalCharge * ChargeScaling; // Apply charge scaling to receptor atoms
                    }
                }

                // Track source residue key/index for diagnostics
                Atom.SourceKey = Pair.Key;
                Atom.SourceIndex = i;

                ReceptorAtoms.Add(Atom);
            }
        }
    }

    // Load ligand atoms (from LigandMap - these are HETATM entries or SDF)
    for (const auto &Pair : ViewerReference->LigandMap)
    {
        const FString &Key = Pair.Key;
        FLigandInfo *LigInfo = Pair.Value;
        if (!LigInfo)
            continue;

        TArray<FMMAtom> Atoms;

        TSet<FString> UnknownElements;
        for (int32 i = 0; i < LigInfo->AtomMeshes.Num(); ++i)
        {
            UStaticMeshComponent *Mesh = LigInfo->AtomMeshes[i];
            if (!Mesh)
                continue;

            FMMAtom Atom;
            Atom.Position = Mesh->GetComponentLocation();
            Atom.bIsReceptor = false;

            // Use element information from LigInfo if available
            if (LigInfo->AtomElements.IsValidIndex(i) && !LigInfo->AtomElements[i].IsEmpty())
                Atom.Element = LigInfo->AtomElements[i].ToUpper();
            else
                Atom.Element = TEXT("C"); // Default

            const FElementData *Data = ElementDatabase.Find(Atom.Element);
            if (Data)
            {
                Atom.Radius = Data->VDWRadius * 100.0f; // Convert to UE units
                Atom.Charge = Data->TypicalCharge;
            }
            else
            {
                if (!UnknownElements.Contains(Atom.Element))
                {
                    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Unknown element '%s' encountered when loading ligand %s. Defaulting to C."), *Atom.Element, *Key);
                    UnknownElements.Add(Atom.Element);
                }
                Atom.Element = TEXT("C");
                const FElementData *Def = ElementDatabase.Find(Atom.Element);
                if (Def)
                {
                    Atom.Radius = Def->VDWRadius * 100.0f;
                    Atom.Charge = Def->TypicalCharge;
                }
            }

            // Track source ligand key/index for diagnostics
            Atom.SourceKey = Key;
            Atom.SourceIndex = i;

            Atoms.Add(Atom);
        }

        AssignPartialCharges(Atoms);
        LigandAtoms.Add(Key, Atoms);
    }
}

FBindingAffinityResult AMMGBSA::CalculateBindingAffinity(const FString &LigandKey)
{
    FBindingAffinityResult Result;
    Result.LigandKey = LigandKey;

    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Request to calculate affinity for %s"), *LigandKey);

    // CHECK FOR PARAMETER CHANGES - ADD THIS LINE
    CheckParameterChanges();

    // Check cache
    if (CachedResults.Contains(LigandKey))
    {
        return CachedResults[LigandKey];
    }

    // Ensure we have atom data loaded
    if (LigandAtoms.Num() == 0 && ViewerReference)
    {
        LoadAtomsFromViewer();
    }

    // Find ligand atoms
    TArray<FMMAtom> *LigandPtr = LigandAtoms.Find(LigandKey);
    if (!LigandPtr || LigandPtr->Num() == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Ligand %s not found or has no atoms"), *LigandKey);
        return Result;
    }

    TArray<FMMAtom> &Ligand = *LigandPtr;
    Result.LigandName = LigandKey;

    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: About to validate structure for %s with %d ligand atoms, %d receptor atoms"),
           *LigandKey, Ligand.Num(), ReceptorAtoms.Num());

    // Validate structure quality before calculation
    if (!ValidateStructureQuality(Ligand, Result))
    {
        UE_LOG(LogTemp, Error, TEXT("MMGBSA: Validation FAILED for %s - bAutoMinimizeOnOverlap=%d"),
               *LigandKey, bAutoMinimizeOnOverlap ? 1 : 0);
        UE_LOG(LogTemp, Error, TEXT("MMGBSA: Quality warning: %s"), *Result.QualityWarning);

        if (bAutoMinimizeOnOverlap)
        {
            UE_LOG(LogTemp, Log, TEXT("MMGBSA: Attempting automatic minimization for %s"), *LigandKey);
            if (MinimizeStructure(LigandKey, 500, 0.1f))
            {
                // Reload pointer after minimization (structure was modified in place)
                LigandPtr = LigandAtoms.Find(LigandKey);
                if (!LigandPtr)
                {
                    UE_LOG(LogTemp, Error, TEXT("MMGBSA: Lost ligand reference after minimization"));
                    Result.bIsValid = false;
                    CachedResults.Add(LigandKey, Result);
                    OnAffinityCalculated.Broadcast(Result);
                    return Result;
                }
                // Update reference to point to newly minimized structure
                TArray<FMMAtom> &LigandMinimized = *LigandPtr;
                // Re-validate after minimization
                if (!ValidateStructureQuality(LigandMinimized, Result))
                {
                    UE_LOG(LogTemp, Error, TEXT("MMGBSA: Structure still has severe overlaps after minimization"));
                    CachedResults.Add(LigandKey, Result);
                    OnAffinityCalculated.Broadcast(Result);
                    return Result;
                }
                UE_LOG(LogTemp, Log, TEXT("MMGBSA: Structure minimized successfully, proceeding with calculation"));
            }
            else
            {
                UE_LOG(LogTemp, Error, TEXT("MMGBSA: Minimization failed for %s"), *LigandKey);
                CachedResults.Add(LigandKey, Result);
                OnAffinityCalculated.Broadcast(Result);
                return Result;
            }
        }
        else
        {
            // Structure is bad and auto-minimize is off - return invalid result
            UE_LOG(LogTemp, Error, TEXT("MMGBSA: Auto-minimize is OFF - returning invalid result for %s"), *LigandKey);
            UE_LOG(LogTemp, Error, TEXT("MMGBSA: To fix: Enable 'Auto Minimize On Overlap' in MMGBSA actor details, or manually minimize the structure"));
            CachedResults.Add(LigandKey, Result);
            OnAffinityCalculated.Broadcast(Result);
            return Result;
        }
    }
    else
    {
        UE_LOG(LogTemp, Log, TEXT("MMGBSA: Structure validation PASSED for %s"), *LigandKey);
    }

    // Create complex (receptor + ligand)
    TArray<FMMAtom> Complex = ReceptorAtoms;
    Complex.Append(Ligand);

    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Calculating for %s (%d atoms)"), *LigandKey, Ligand.Num());

    // Diagnostic: list ligand atoms with element/position/charge/radius
    float TotalAbsCharge = 0.0f;
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Atoms for %s:"), *LigandKey);
    for (int32 i = 0; i < Ligand.Num(); ++i)
    {
        const FMMAtom &A = Ligand[i];
        UE_LOG(LogTemp, Log, TEXT("  Atom %d: Element=%s Pos=(%.2f, %.2f, %.2f) Charge=%.4f Radius=%.2f"), i, *A.Element, A.Position.X, A.Position.Y, A.Position.Z, A.Charge, A.Radius);
        TotalAbsCharge += FMath::Abs(A.Charge);
    }
    if (FMath::IsNearlyZero(TotalAbsCharge, 1e-6f))
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Total absolute ligand charge for %s is ~0.0; charges may not have been assigned correctly."), *LigandKey);
    }

    // Per-atom VDW overlap checks (ligand vs receptor) with detailed info
    TArray<FMMOverlapInfo> OverlapList;
    const int32 MaxLogPairs = 10;
    int32 LoggedPairs = 0;
    const float SevereFraction = 0.5f; // Dist < 0.5*(ri + rj) => severe overlap
    for (int32 i = 0; i < ReceptorAtoms.Num(); ++i)
    {
        const FMMAtom &R = ReceptorAtoms[i];
        for (int32 j = 0; j < Ligand.Num(); ++j)
        {
            const FMMAtom &L = Ligand[j];
            float DistUE = FVector::Dist(R.Position, L.Position);
            float DistA = DistUE / 100.0f; // Å
            float RiA = R.Radius / 100.0f;
            float RjA = L.Radius / 100.0f;
            float SumR = RiA + RjA;

            if (DistA < SevereFraction * SumR)
            {
                FMMOverlapInfo O;
                O.ReceptorSourceKey = R.SourceKey;
                O.ReceptorAtomIndex = R.SourceIndex;
                O.ReceptorElement = R.Element;
                O.ReceptorPosition = R.Position;
                O.LigandAtomIndex = L.SourceIndex;
                O.LigandElement = L.Element;
                O.LigandPosition = L.Position;
                O.DistanceA = DistA;
                O.SumRadiiA = SumR;

                // Estimate pair VDW energy for context
                float Epsilon = FMath::Sqrt(GetVDWEpsilon(R.Element) * GetVDWEpsilon(L.Element));
                float Sigma = (RiA + RjA) * 0.5f;
                float DistUse = DistA;
                if (bUseSoftCoreVDW)
                {
                    DistUse = FMath::Sqrt(DistA * DistA + SoftCoreAlpha * SoftCoreAlpha);
                }
                float SigmaOverR = Sigma / FMath::Max(0.0001f, DistUse);
                float SigmaOverR6 = FMath::Pow(SigmaOverR, 6.0f);
                float SigmaOverR12 = SigmaOverR6 * SigmaOverR6;
                float PairE = 4.0f * Epsilon * (SigmaOverR12 - SigmaOverR6);
                if (!FMath::IsFinite(PairE))
                    PairE = FMath::Sign(PairE) * MaxPairVDW;
                PairE = FMath::Clamp(PairE, -MaxPairVDW, MaxPairVDW);
                O.PairVDW = PairE;

                OverlapList.Add(O);

                if (LoggedPairs < MaxLogPairs)
                {
                    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Severe overlap between receptor %s atom %d (%s) and ligand atom %d (%s): Dist=%.3f Å SumR=%.3f Å PairVDW=%.3f kcal/mol"),
                           *O.ReceptorSourceKey, O.ReceptorAtomIndex, *O.ReceptorElement, O.LigandAtomIndex, *O.LigandElement, O.DistanceA, O.SumRadiiA, O.PairVDW);
                    LoggedPairs++;
                }
            }
        }
    }

    Result.OverlapCount = OverlapList.Num();
    Result.bHasSevereOverlap = (Result.OverlapCount > 0);
    Result.Overlaps = OverlapList;
    if (Result.bHasSevereOverlap)
    {
        if (ViewerReference)
        {
            ViewerReference->HighlightOverlapAtoms(Result.Overlaps);
        }

        if (!bUseSoftCoreVDW)
        {
            UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Found %d severe ligand/receptor VDW overlaps for %s; aborting calculation and marking result invalid."), Result.OverlapCount, *LigandKey);

            // Mark as invalid and return early so UI shows a failure instead of misleading energies
            Result.bIsValid = false;
            CachedResults.Add(LigandKey, Result);
            OnAffinityCalculated.Broadcast(Result);
            return Result;
        }
        else
        {
            UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Found %d severe ligand/receptor VDW overlaps for %s; soft-core VDW enabled, proceeding with softened pair potentials."), Result.OverlapCount, *LigandKey);
            // Continue with calculation but keep bHasSevereOverlap true so UI can surface it
        }
    }

    // Calculate MM energies using ONLY pairwise receptor-ligand interactions
    // This avoids catastrophic cancellation from large internal energies
    Result.DeltaEMM = CalculateMMEnergyPair(ReceptorAtoms, Ligand);

    // Still calculate full energies for diagnostic logging, but don't use them
    Result.EMM_Complex = CalculateMMEnergy(Complex);
    Result.EMM_Receptor = CalculateMMEnergy(ReceptorAtoms);
    Result.EMM_Ligand = CalculateMMEnergy(Ligand);

    UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: ΔEMM (receptor-ligand pairs only) = %.6f kcal/mol"), Result.DeltaEMM);

    // Warn if internal energies are unrealistically large (indicates structural problems)
    if (FMath::Abs(Result.EMM_Complex) > 50000.0f || FMath::Abs(Result.EMM_Receptor) > 50000.0f)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Very large internal MM energies detected (Complex=%.0f, Receptor=%.0f) - structure may need energy minimization"),
               Result.EMM_Complex, Result.EMM_Receptor);
    }

    // Calculate GB solvation energies
    // Use a more numerically stable approach: compute receptor-ligand GB interactions
    // and apply empirical scaling to account for desolvation and other effects
    TArray<FMMAtom> ReceptorCopy = ReceptorAtoms;
    TArray<FMMAtom> LigandCopy = Ligand;
    TArray<FMMAtom> ComplexCopy = Complex;

    // Calculate SASA early so we can use it for desolvation estimate
    Result.SA_Complex = CalculateSurfaceAreaEnergy(ComplexCopy);
    Result.SA_Receptor = CalculateSurfaceAreaEnergy(ReceptorCopy);
    Result.SA_Ligand = CalculateSurfaceAreaEnergy(LigandCopy);
    Result.DeltaGSA = Result.SA_Complex - Result.SA_Receptor - Result.SA_Ligand;

    // Calculate full GB energies for logging/debugging
    Result.GB_Complex = CalculateGBEnergy(ComplexCopy);
    Result.GB_Receptor = CalculateGBEnergy(ReceptorCopy);
    Result.GB_Ligand = CalculateGBEnergy(LigandCopy);

    // Compute ΔGGB using pairwise receptor-ligand interactions
    // This avoids the large cancellation from full complex/receptor/ligand calculation
    float PairwiseGB = CalculateGBEnergyPair(ReceptorAtoms, Ligand);

    UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: Raw pairwise GB = %.6f kcal/mol"), PairwiseGB);

    // Estimate desolvation penalty based on buried surface area
    // When receptor and ligand bind, they lose solvation, which is unfavorable
    // This partially offsets the favorable receptor-ligand GB interactions
    // DeltaSASA is negative (buried surface), so we need to convert to positive penalty
    // Desolvation penalty for electrostatic solvation is smaller than for nonpolar
    // Use a conservative estimate: ~0.005 kcal/mol per Å² of buried surface
    float BuriedSASA = -(Result.SA_Complex - Result.SA_Receptor - Result.SA_Ligand); // Make positive
    float DesolvationPenalty = BuriedSASA * 0.005f;                                  // Positive penalty (conservative estimate)

    // Apply empirical scaling to pairwise GB interactions
    // The scaling accounts for approximations in the GB model and charge assignments
    float ScaledPairwiseGB = PairwiseGB * GBScale;

    // Net GB energy change = scaled pairwise interactions - desolvation penalty
    Result.DeltaGGB = ScaledPairwiseGB + DesolvationPenalty;

    UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: Pairwise GB=%.3f, scaled=%.3f, desolvation=%.3f, net ΔGGB=%.3f"),
           PairwiseGB, ScaledPairwiseGB, DesolvationPenalty, Result.DeltaGGB);

    if (!FMath::IsNearlyEqual(GBScale, 1.0f))
    {
        UE_LOG(LogTemp, Log, TEXT("MMGBSA: Applied GBScale=%.3f to ΔGGB (raw=%.3f -> scaled=%.3f, with desolvation=%.3f)"),
               GBScale, PairwiseGB, ScaledPairwiseGB, Result.DeltaGGB);
    }

    // Clamp ΔGGB magnitude to protect against numerical extremes
    if (!FMath::IsFinite(Result.DeltaGGB))
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: ΔGGB computed non-finite (NaN/inf); clamping to sign-preserving finite value."));
        Result.DeltaGGB = FMath::Sign(Result.DeltaGGB) * MaxDeltaGB;
    }
    if (FMath::Abs(Result.DeltaGGB) > MaxDeltaGB)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: ΔGGB magnitude (%.3f) exceeds MaxDeltaGB (%.3f); clamping."), Result.DeltaGGB, MaxDeltaGB);
        Result.DeltaGGB = FMath::Clamp(Result.DeltaGGB, -MaxDeltaGB, MaxDeltaGB);
    }

    // Surface area terms already calculated above for desolvation estimate

    // Total binding free energy
    Result.DeltaG_Binding = Result.DeltaEMM + Result.DeltaGGB + Result.DeltaGSA;

    // Sanity check: if ΔEMM is unreasonably large (> 100 kcal/mol), it's likely due to severe overlaps
    // In these cases, the calculation is unreliable and should be flagged
    if (FMath::Abs(Result.DeltaEMM) > 100.0f && !Result.bHasSevereOverlap)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: ΔEMM magnitude (%.1f) is very large but no severe overlaps detected - possible numerical instability"), Result.DeltaEMM);
    }

    // Clamp final ΔG to reasonable range (-50 to +50 kcal/mol)
    // Values outside this range are almost certainly numerical artifacts
    const float MaxDeltaG = 50.0f;
    if (!FMath::IsFinite(Result.DeltaG_Binding))
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: ΔG computed as non-finite; marking result invalid"));
        Result.DeltaG_Binding = 0.0f;
        Result.bIsValid = false;
    }
    else if (FMath::Abs(Result.DeltaG_Binding) > MaxDeltaG)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: ΔG magnitude (%.1f) exceeds reasonable range (±%.0f kcal/mol); clamping and marking uncertain"),
               Result.DeltaG_Binding, MaxDeltaG);
        Result.DeltaG_Binding = FMath::Clamp(Result.DeltaG_Binding, -MaxDeltaG, MaxDeltaG);
        // Don't mark invalid, but user should be aware this is clamped
    }

    // Convert to Ki
    Result.Ki_uM = DeltaGToKi(Result.DeltaG_Binding);
    Result.AffinityClass = ClassifyAffinity(Result.DeltaG_Binding);

    // Log breakdown of components for diagnostics
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: %s | ΔG = %.2f kcal/mol | Ki = %.2f µM | %s"), *LigandKey, Result.DeltaG_Binding, Result.Ki_uM, *Result.AffinityClass);
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Energy breakdown for %s:"), *LigandKey);
    UE_LOG(LogTemp, Log, TEXT("  ΔEMM (receptor-ligand) = %.3f kcal/mol"), Result.DeltaEMM);
    UE_LOG(LogTemp, Log, TEXT("  ΔGGB (solvation)       = %.3f kcal/mol"), Result.DeltaGGB);
    UE_LOG(LogTemp, Log, TEXT("  ΔGSA (surface area)    = %.3f kcal/mol"), Result.DeltaGSA);
    UE_LOG(LogTemp, Verbose, TEXT("  Internal energies (diagnostic): EMM_Complex=%.1f EMM_Receptor=%.1f EMM_Ligand=%.1f"),
           Result.EMM_Complex, Result.EMM_Receptor, Result.EMM_Ligand);

    // Warn on unusually large energy components which likely indicate numeric instability or atom overlap
    const float WarnThreshold = 1000.0f; // kcal/mol
    if (FMath::Abs(Result.DeltaEMM) > WarnThreshold || FMath::Abs(Result.DeltaGGB) > WarnThreshold || FMath::Abs(Result.DeltaGSA) > WarnThreshold)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Large energy component detected for %s (ΔEMM=%.3f ΔGB=%.3f ΔSA=%.3f); check for overlapping atoms or incorrect radii/charges"), *LigandKey, Result.DeltaEMM, Result.DeltaGGB, Result.DeltaGSA);
    }

    // If ΔG is (unexpectedly) zero but atoms were present, warn about possible missing element/charge data
    if (Ligand.Num() > 0 && FMath::IsNearlyZero(Result.DeltaG_Binding, 1e-4f))
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: ΔG for %s computed as ~0.0 with components EMM=%.3f GB=%.3f SA=%.3f; verify atom elements and charges are assigned correctly"),
               *LigandKey, Result.DeltaEMM, Result.DeltaGGB, Result.DeltaGSA);
    }

    // Mark valid and cache result
    Result.bIsValid = true;
    CachedResults.Add(LigandKey, Result);

    // Broadcast event
    OnAffinityCalculated.Broadcast(Result);

    return Result;
}

TArray<FBindingAffinityResult> AMMGBSA::CalculateAllBindingAffinities()
{
    TArray<FBindingAffinityResult> Results;

    for (const auto &Pair : LigandAtoms)
    {
        FBindingAffinityResult Result = CalculateBindingAffinity(Pair.Key);
        Results.Add(Result);
    }

    // Sort by binding affinity (most negative = strongest)
    Results.Sort([](const FBindingAffinityResult &A, const FBindingAffinityResult &B)
                 { return A.DeltaG_Binding < B.DeltaG_Binding; });

    return Results;
}

bool AMMGBSA::GetCachedResult(const FString &LigandKey, FBindingAffinityResult &OutResult) const
{
    const FBindingAffinityResult *Found = CachedResults.Find(LigandKey);
    if (Found)
    {
        OutResult = *Found;
        return true;
    }
    return false;
}

// MM Energy = Electrostatic + Van der Waals
float AMMGBSA::CalculateMMEnergy(const TArray<FMMAtom> &Atoms) const
{
    float Elec = CalculateElectrostaticEnergy(Atoms);
    float VDW = CalculateVanDerWaalsEnergy(Atoms);
    return Elec + VDW;
}

// Pairwise MM energy between two groups (used for receptor-ligand interactions to avoid large cancellation)
float AMMGBSA::CalculateMMEnergyPair(const TArray<FMMAtom> &Atoms1, const TArray<FMMAtom> &Atoms2) const
{
    float Energy = 0.0f;
    const float CutoffA = this->ReceptorLigandCutoffA;

    for (int32 i = 0; i < Atoms1.Num(); ++i)
    {
        const FMMAtom &A = Atoms1[i];
        for (int32 j = 0; j < Atoms2.Num(); ++j)
        {
            const FMMAtom &B = Atoms2[j];
            float Dist = FVector::Dist(A.Position, B.Position) / 100.0f; // Å

            // MODIFIED: More aggressive minimum distance clamping
            if (Dist < 1.5f) // INCREASED from 1.0f
            {
                Dist = 1.5f; // Clamp to more realistic minimum
                static int32 MinDistWarnings = 0;
                if (MinDistWarnings < 3)
                {
                    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Severe atomic overlap (dist < 1.5 Å), clamping to minimum"));
                    MinDistWarnings++;
                }
            }

            if (Dist > CutoffA)
                continue;

            // Distance-dependent dielectric
            float PairDielectric = FMath::Max(4.0f, 4.0f * Dist);
            float CoulombKPair = 332.0636f / PairDielectric;
            float Elec = CoulombKPair * A.Charge * B.Charge / Dist;

            // VDW with soft-core
            float Epsilon = FMath::Sqrt(GetVDWEpsilon(A.Element) * GetVDWEpsilon(B.Element));
            float Sigma = (GetVDWRadius(A.Element) + GetVDWRadius(B.Element)) * 0.5f;

            float DistUse = Dist;
            if (bUseSoftCoreVDW)
            {
                DistUse = FMath::Sqrt(Dist * Dist + SoftCoreAlpha * SoftCoreAlpha);
            }

            float SigmaOverR = Sigma / FMath::Max(0.0001f, DistUse);
            float SigmaOverR6 = FMath::Pow(SigmaOverR, 6.0f);
            float SigmaOverR12 = SigmaOverR6 * SigmaOverR6;
            float VDW = 4.0f * Epsilon * (SigmaOverR12 - SigmaOverR6);
            if (!FMath::IsFinite(VDW))
                VDW = FMath::Sign(VDW) * MaxPairVDW;
            VDW = FMath::Clamp(VDW, -MaxPairVDW, MaxPairVDW);

            // MODIFIED: Stronger distance-dependent scaling
            float DistScale = 1.0f;
            if (Dist < 3.0f)
            {
                DistScale = FMath::Lerp(0.05f, 1.0f, FMath::Clamp((Dist - 0.5f) / 2.5f, 0.0f, 1.0f));
            }
            else if (Dist > 5.0f)
            {
                DistScale = FMath::Lerp(1.0f, 0.6f, FMath::Clamp((Dist - 5.0f) / (CutoffA - 5.0f), 0.0f, 1.0f));
            }

            Energy += (Elec + VDW) * DistScale;
        }
    }

    return Energy;
}

float AMMGBSA::CalculateElectrostaticEnergy(const TArray<FMMAtom> &Atoms) const
{
    // Coulomb's law: E = (1/4πε₀) * Σᵢⱼ (qᵢqⱼ / rᵢⱼ)
    // In vacuum: k = 332.0636 kcal·Å/(mol·e²)
    const float CoulombK = 332.0636f / InteriorDielectric;
    float Energy = 0.0f;
    const float CutoffA = 12.0f; // Å - ignore interactions beyond this distance

    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        for (int32 j = i + 1; j < Atoms.Num(); ++j)
        {
            float Dist = FVector::Dist(Atoms[i].Position, Atoms[j].Position) / 100.0f; // Convert to Å
            if (Dist < 0.5f)
                continue; // Avoid singularity
            if (Dist > CutoffA)
                continue; // Skip long-range pairs to avoid huge cumulative sums

            float E = CoulombK * Atoms[i].Charge * Atoms[j].Charge / Dist;
            Energy += E;
        }
    }

    return Energy;
}

float AMMGBSA::CalculateVanDerWaalsEnergy(const TArray<FMMAtom> &Atoms) const
{
    // Lennard-Jones: E = Σᵢⱼ 4εᵢⱼ[(σᵢⱼ/rᵢⱼ)¹² - (σᵢⱼ/rᵢⱼ)⁶]
    float Energy = 0.0f;
    const float CutoffA = 12.0f; // Å - ignore interactions beyond this distance

    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        for (int32 j = i + 1; j < Atoms.Num(); ++j)
        {
            float Dist = FVector::Dist(Atoms[i].Position, Atoms[j].Position) / 100.0f; // Convert to Å
            if (Dist < 0.5f)
                continue;
            if (Dist > CutoffA)
                continue; // Skip long-range VDW pairs

            // Lorentz-Berthelot combining rules
            float Epsilon = FMath::Sqrt(
                GetVDWEpsilon(Atoms[i].Element) *
                GetVDWEpsilon(Atoms[j].Element));

            float Sigma = (GetVDWRadius(Atoms[i].Element) +
                           GetVDWRadius(Atoms[j].Element)) *
                          0.5f;

            float DistUse = Dist;
            if (bUseSoftCoreVDW)
            {
                DistUse = FMath::Sqrt(Dist * Dist + SoftCoreAlpha * SoftCoreAlpha);
            }

            float SigmaOverR = Sigma / FMath::Max(0.0001f, DistUse);
            float SigmaOverR6 = FMath::Pow(SigmaOverR, 6.0f);
            float SigmaOverR12 = SigmaOverR6 * SigmaOverR6;

            float E = 4.0f * Epsilon * (SigmaOverR12 - SigmaOverR6);
            if (!FMath::IsFinite(E))
                E = FMath::Sign(E) * MaxPairVDW;
            E = FMath::Clamp(E, -MaxPairVDW, MaxPairVDW);
            Energy += E;
        }
    }

    return Energy;
}

// Generalized Born implicit solvation
float AMMGBSA::CalculateGBEnergy(const TArray<FMMAtom> &Atoms)
{
    TArray<FMMAtom> AtomsCopy = Atoms;
    CalculateBornRadii(AtomsCopy);
    const float Factor = -0.5f * 332.0636f * (1.0f / this->InteriorDielectric - 1.0f / this->ExteriorDielectric);
    float Energy = 0.0f;

    // Per-atom self-energy with clamping to avoid huge spikes
    for (const FMMAtom &Atom : AtomsCopy)
    {
        if (Atom.GBRadius > 0.0f)
        {
            float PerAtom = Factor * Atom.Charge * Atom.Charge / Atom.GBRadius;
            if (!FMath::IsFinite(PerAtom))
                PerAtom = FMath::Sign(PerAtom) * MaxGBPerAtom;
            float Clamped = FMath::Clamp(PerAtom, -MaxGBPerAtom, MaxGBPerAtom);
            if (!FMath::IsNearlyEqual(PerAtom, Clamped))
            {
                UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: Clamped per-atom GB contribution from %.3f to %.3f for element %s"), PerAtom, Clamped, *Atom.Element);
            }
            Energy += Clamped;
        }
    }

    const float CutoffA = this->ReceptorLigandGBCutoffA;
    for (int32 i = 0; i < AtomsCopy.Num(); ++i)
    {
        for (int32 j = i + 1; j < AtomsCopy.Num(); ++j)
        {
            float Dist = FVector::Dist(AtomsCopy[i].Position, AtomsCopy[j].Position) / 100.0f;

            // Only include pairwise GB contribution if within cutoff
            if (Dist > CutoffA)
                continue;

            float Ai = AtomsCopy[i].GBRadius;
            float Aj = AtomsCopy[j].GBRadius;
            if (Ai < 0.5f || Aj < 0.5f)
            {
                static int32 RadiiWarnings = 0;
                if (RadiiWarnings < 3)
                {
                    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Unrealistic Born radius: Ai=%.3f Aj=%.3f, skipping pair"), Ai, Aj);
                    RadiiWarnings++;
                }
                continue;
            }
            float AiAj = Ai * Aj;
            float Exp = FMath::Exp(-Dist * Dist / (4.0f * AiAj));
            float fGB_arg = Dist * Dist + AiAj * Exp;
            if (fGB_arg < 0.01f)
            {
                fGB_arg = 0.01f;
            }
            float fGB = 1.0f / FMath::Sqrt(fGB_arg);
            fGB = FMath::Min(fGB, 100.0f);

            float PairTerm = 2.0f * Factor * AtomsCopy[i].Charge * AtomsCopy[j].Charge * fGB;
            if (!FMath::IsFinite(PairTerm))
                PairTerm = FMath::Sign(PairTerm) * MaxGBPair;
            PairTerm = FMath::Clamp(PairTerm, -MaxGBPair, MaxGBPair);
            if (!FMath::IsNearlyZero(PairTerm) && FMath::Abs(PairTerm) == MaxGBPair)
            {
                UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: Clamped GB pair term between atoms %d and %d to %.3f"), i, j, PairTerm);
            }

            Energy += PairTerm;
        }
    }

    return Energy;
}

// Pairwise GB energy between two groups (used for receptor-ligand interactions to avoid large cancellation)
// This computes the GB energy change due to receptor-ligand interactions, avoiding the numerical
// instability from computing full GB_Complex - GB_Receptor - GB_Ligand with large self-energy terms
float AMMGBSA::CalculateGBEnergyPair(const TArray<FMMAtom> &Atoms1, const TArray<FMMAtom> &Atoms2)
{
    TArray<FMMAtom> Combined = Atoms1;
    Combined.Append(Atoms2);
    CalculateBornRadii(Combined);

    TArray<FMMAtom> Atoms1Copy = Atoms1;
    TArray<FMMAtom> Atoms2Copy = Atoms2;
    for (int32 i = 0; i < Atoms1Copy.Num(); ++i)
    {
        Atoms1Copy[i].GBRadius = Combined[i].GBRadius;
    }
    for (int32 i = 0; i < Atoms2Copy.Num(); ++i)
    {
        Atoms2Copy[i].GBRadius = Combined[Atoms1.Num() + i].GBRadius;
    }

    const float Factor = -0.5f * 332.0636f * (1.0f / InteriorDielectric - 1.0f / ExteriorDielectric);
    float Energy = 0.0f;
    const float CutoffA = ReceptorLigandGBCutoffA;

    const float DistScaleStart = 3.0f;
    const float DistScaleEnd = CutoffA;

    for (int32 i = 0; i < Atoms1Copy.Num(); ++i)
    {
        const FMMAtom &A = Atoms1Copy[i];
        for (int32 j = 0; j < Atoms2Copy.Num(); ++j)
        {
            const FMMAtom &B = Atoms2Copy[j];
            float Dist = FVector::Dist(A.Position, B.Position) / 100.0f;

            if (Dist > CutoffA || Dist < 0.5f)
                continue;

            float Ai = A.GBRadius;
            float Aj = B.GBRadius;

            // NEW: Skip unrealistic Born radii
            if (Ai < 0.5f || Aj < 0.5f)
            {
                static int32 RadiiWarnings = 0;
                if (RadiiWarnings < 3)
                {
                    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Unrealistic Born radius: Ai=%.3f Aj=%.3f, skipping pair"), Ai, Aj);
                    RadiiWarnings++;
                }
                continue;
            }

            float AiAj = Ai * Aj;
            float Exp = FMath::Exp(-Dist * Dist / (4.0f * AiAj));

            // NEW: Prevent division by near-zero
            float fGB_arg = Dist * Dist + AiAj * Exp;
            if (fGB_arg < 0.01f)
            {
                fGB_arg = 0.01f;
            }
            float fGB = 1.0f / FMath::Sqrt(fGB_arg);

            // NEW: Clamp fGB itself to prevent extreme values
            fGB = FMath::Min(fGB, 100.0f);

            float PairTerm = 2.0f * Factor * A.Charge * B.Charge * fGB;

            // Distance-dependent scaling
            float DistScale = 1.0f;
            if (Dist > DistScaleStart)
            {
                float T = (Dist - DistScaleStart) / (DistScaleEnd - DistScaleStart);
                DistScale = FMath::Lerp(1.0f, 0.3f, FMath::Clamp(T, 0.0f, 1.0f));
            }
            PairTerm *= DistScale;

            if (!FMath::IsFinite(PairTerm))
                PairTerm = FMath::Sign(PairTerm) * MaxGBPair;
            PairTerm = FMath::Clamp(PairTerm, -MaxGBPair, MaxGBPair);

            Energy += PairTerm;
        }
    }

    return Energy;
}

void AMMGBSA::CalculateBornRadii(TArray<FMMAtom> &Atoms)
{
    // Simplified Born radii calculation
    // More accurate methods: Still, Onufriev-Bashford-Case, etc.

    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        float Rho = GetEffectiveBornRadius(Atoms[i], Atoms);

        const FElementData *Data = ElementDatabase.Find(Atoms[i].Element);
        float Scale = Data ? Data->GBRadiusScale : 1.0f;

        Atoms[i].GBRadius = Atoms[i].Radius * Scale / 100.0f; // Convert to Å

        // Adjust based on burial
        if (Rho > 0.0f)
        {
            Atoms[i].GBRadius *= (1.0f + 0.5f * FMath::Exp(-Rho / 5.0f));
        }

        // Enforce minimum GB radius to avoid tiny denominators
        if (Atoms[i].GBRadius < MinGBRadius)
        {
            UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: Bumped GB radius for element %s from %.3f to minimum %.3f Å"), *Atoms[i].Element, Atoms[i].GBRadius, MinGBRadius);
            Atoms[i].GBRadius = MinGBRadius;
        }
    }
}

float AMMGBSA::GetEffectiveBornRadius(const FMMAtom &Atom, const TArray<FMMAtom> &AllAtoms)
{
    // Count nearby atoms (burial factor)
    float Rho = 0.0f;
    const float CutoffSq = FMath::Square(1000.0f); // 10 Å

    for (const FMMAtom &Other : AllAtoms)
    {
        if (&Other == &Atom)
            continue;

        float DistSq = FVector::DistSquared(Atom.Position, Other.Position);
        if (DistSq < CutoffSq && DistSq > 1.0f)
        {
            float Dist = FMath::Sqrt(DistSq) / 100.0f; // Å
            Rho += FMath::Exp(-Dist / 3.0f);           // Exponential decay
        }
    }

    return Rho;
}

float AMMGBSA::CalculateSurfaceAreaEnergy(const TArray<FMMAtom> &Atoms)
{
    TArray<FMMAtom> AtomsCopy = Atoms;
    CalculateSASA(AtomsCopy);
    float TotalSASA = 0.0f;
    for (const FMMAtom &Atom : AtomsCopy)
    {
        TotalSASA += Atom.SASA;
    }
    return this->SurfaceTension * TotalSASA;
}

void AMMGBSA::CalculateSASA(TArray<FMMAtom> &Atoms)
{
    // Simplified SASA calculation using sphere overlap
    // More accurate: Lee-Richards algorithm, Shrake-Rupley

    const int32 NumSpherePoints = 92;                   // Fibonacci sphere sampling
    const float ProbeRad = SolventProbeRadius * 100.0f; // Convert to UE units

    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        FMMAtom &Atom = Atoms[i];
        float ExtendedRadius = Atom.Radius + ProbeRad;
        int32 AccessiblePoints = 0;

        // Generate test points on sphere surface
        for (int32 p = 0; p < NumSpherePoints; ++p)
        {
            // Fibonacci sphere
            float Phi = FMath::Acos(1.0f - 2.0f * (p + 0.5f) / NumSpherePoints);
            float Theta = PI * (1.0f + FMath::Sqrt(5.0f)) * p;

            FVector TestPoint = Atom.Position + ExtendedRadius * FVector(
                                                                     FMath::Sin(Phi) * FMath::Cos(Theta),
                                                                     FMath::Sin(Phi) * FMath::Sin(Theta),
                                                                     FMath::Cos(Phi));

            // Check if point is accessible (not inside another atom)
            bool bAccessible = true;
            for (int32 j = 0; j < Atoms.Num(); ++j)
            {
                if (i == j)
                    continue;

                float Dist = FVector::Dist(TestPoint, Atoms[j].Position);
                if (Dist < (Atoms[j].Radius + ProbeRad))
                {
                    bAccessible = false;
                    break;
                }
            }

            if (bAccessible)
                AccessiblePoints++;
        }

        // SASA = 4πr² * (accessible fraction), convert to Å³
        float Fraction = (float)AccessiblePoints / NumSpherePoints;
        Atom.SASA = 4.0f * PI * ExtendedRadius * ExtendedRadius * Fraction / 10000.0f; // To Å³
    }
}

void AMMGBSA::AssignPartialCharges(TArray<FMMAtom> &Atoms)
{
    // Improved charge assignment based on element type and chemical context
    // This is still a simplified model - real charges would come from QM or force field

    TSet<FString> WarnedElements;

    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        FMMAtom &Atom = Atoms[i];
        float BaseCharge = 0.0f;

        // Assign charges based on element with chemical context
        FString El = Atom.Element.ToUpper();

        if (El == TEXT("C"))
        {
            // Carbons: small positive charge in polar environments, neutral otherwise
            int32 NearbyHeteroatoms = 0;
            for (int32 j = 0; j < Atoms.Num(); ++j)
            {
                if (i == j)
                    continue;
                FString OtherEl = Atoms[j].Element.ToUpper();
                if (OtherEl == TEXT("N") || OtherEl == TEXT("O") || OtherEl == TEXT("S"))
                {
                    float Dist = FVector::Dist(Atom.Position, Atoms[j].Position) / 100.0f; // Å
                    if (Dist < 2.5f)
                    {
                        NearbyHeteroatoms++;
                    }
                }
            }

            if (NearbyHeteroatoms > 0)
            {
                BaseCharge = 0.0375f * NearbyHeteroatoms;  // REDUCED from 0.075
                BaseCharge = FMath::Min(BaseCharge, 0.1f); // REDUCED cap from 0.2
            }
            else
            {
                BaseCharge = 0.0125f; // REDUCED from 0.025
            }
        }
        else if (El == TEXT("N"))
        {
            BaseCharge = -0.1f; // REDUCED from -0.2
        }
        else if (El == TEXT("O"))
        {
            BaseCharge = -0.125f; // REDUCED from -0.25
        }
        else if (El == TEXT("S"))
        {
            BaseCharge = -0.05f; // REDUCED from -0.1
        }
        else if (El == TEXT("P"))
        {
            BaseCharge = 0.2f; // REDUCED from 0.4
        }
        else if (El == TEXT("H"))
        {
            BaseCharge = 0.025f; // REDUCED from 0.05
        }
        else if (El == TEXT("F") || El == TEXT("CL") || El == TEXT("BR") || El == TEXT("I"))
        {
            BaseCharge = -0.05f; // REDUCED from -0.1
        }
        else
        {
            // Metals and others: use database value or default to neutral
            const FElementData *Data = ElementDatabase.Find(El);
            BaseCharge = Data ? Data->TypicalCharge : 0.0f;
        }

        // Apply ChargeScaling to modulate electrostatic strength
        Atom.Charge = BaseCharge * ChargeScaling;

        // Warn about zero charges only for unexpected cases
        if (FMath::IsNearlyZero(Atom.Charge, 1e-6f) && !WarnedElements.Contains(Atom.Element))
        {
            if (ChargeScaling < 0.001f)
            {
                // If ChargeScaling is essentially zero, this is expected
                UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: Charges scaled to ~0 due to ChargeScaling=%.4f"), ChargeScaling);
            }
            else if (El != TEXT("C") && El != TEXT("H"))
            {
                // Warn if heteroatoms have zero charge unexpectedly
                UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Element '%s' has zero charge after assignment (BaseCharge=%.3f, ChargeScaling=%.3f)"),
                       *Atom.Element, BaseCharge, ChargeScaling);
            }
            WarnedElements.Add(Atom.Element);
        }
    }

    // Neutralize the molecule by distributing excess charge
    // This prevents unphysical long-range electrostatics
    // Only neutralize if total charge is unreasonably large (> 2.0 e)
    // Small net charges are physically reasonable and should be preserved
    float TotalCharge = 0.0f;
    for (const FMMAtom &A : Atoms)
    {
        TotalCharge += A.Charge;
    }

    const float MaxAllowedCharge = 2.0f; // Only neutralize beyond this threshold

    if (FMath::Abs(TotalCharge) > MaxAllowedCharge)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Total charge %.4f e exceeds threshold %.1f e, neutralizing"),
               TotalCharge, MaxAllowedCharge);

        // Distribute excess charge evenly across all atoms
        float Correction = -TotalCharge / Atoms.Num();
        for (FMMAtom &A : Atoms)
        {
            A.Charge += Correction;
        }

        // Verify neutralization
        float FinalCharge = 0.0f;
        for (const FMMAtom &A : Atoms)
        {
            FinalCharge += A.Charge;
        }
        UE_LOG(LogTemp, Log, TEXT("MMGBSA: Total charge after neutralization: %.4f e"), FinalCharge);
    }
    else
    {
        UE_LOG(LogTemp, Log, TEXT("MMGBSA: Total charge %.4f e is acceptable (< %.1f e), preserving charges"),
               TotalCharge, MaxAllowedCharge);
    }
}

float AMMGBSA::EstimateChargeFromElement(const FString &Element) const
{
    const FElementData *Data = ElementDatabase.Find(Element.ToUpper());
    return Data ? Data->TypicalCharge : 0.0f;
}

float AMMGBSA::GetVDWRadius(const FString &Element) const
{
    const FElementData *Data = ElementDatabase.Find(Element.ToUpper());
    return Data ? Data->VDWRadius : 1.70f;
}

float AMMGBSA::GetVDWEpsilon(const FString &Element) const
{
    const FElementData *Data = ElementDatabase.Find(Element.ToUpper());
    return Data ? Data->VDWEpsilon : 0.086f;
}

void AMMGBSA::OnViewerLigandsLoaded()
{
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Viewer reported ligands loaded, refreshing atom cache"));
    LoadAtomsFromViewer();
}

float AMMGBSA::DeltaGToKi(float DeltaG_kcal_mol) const
{
    // ΔG = RT ln(Ki)  => Ki = exp(ΔG / RT)
    const double R = 0.001987;
    const double RT = R * (double)Temperature;
    if (RT <= 0.0)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Invalid temperature RT <= 0 in DeltaGToKi"));
        return FLT_MAX;
    }

    const double Exponent = (double)DeltaG_kcal_mol / RT;

    // Guard against overflow/underflow when exponent is extreme
    double Ki_M = 0.0;
    bool bOverflow = false;
    if (Exponent > 700.0) // exp(700) is already ~1e304, near double overflow
    {
        bOverflow = true;
        Ki_M = DBL_MAX;
    }
    else if (Exponent < -700.0)
    {
        // Extremely favorable binding -> Ki effectively zero
        Ki_M = 0.0;
    }
    else
    {
        Ki_M = FMath::Exp((float)Exponent);
        if (!FMath::IsFinite((float)Ki_M))
            bOverflow = true;
    }

    double Ki_uM = Ki_M * 1.0e6;
    if (!FMath::IsFinite((float)Ki_uM) || bOverflow || Ki_uM > 1e20)
    {
        // Clamp to a large finite value for UI/display rather than returning inf
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Ki overflow/clamped for ΔG=%.3f kcal/mol (exponent=%.3f). Returning large finite Ki."), DeltaG_kcal_mol, Exponent);
        return 1e20f; // 1e20 µM ~ effectively infinite
    }

    return (float)Ki_uM;
}

FString AMMGBSA::ClassifyAffinity(float DeltaG_kcal_mol) const
{
    // Classification based on ΔG
    // ΔG < -10: Very strong (sub-nM)
    // -10 to -7: Strong (nM to low µM)
    // -7 to -5: Moderate (µM)
    // -5 to -3: Weak (high µM to mM)
    // > -3: Very weak (mM+)

    if (DeltaG_kcal_mol < -10.0f)
        return TEXT("Very Strong (pM-nM)");
    else if (DeltaG_kcal_mol < -7.0f)
        return TEXT("Strong (nM-µM)");
    else if (DeltaG_kcal_mol < -5.0f)
        return TEXT("Moderate (µM)");
    else if (DeltaG_kcal_mol < -3.0f)
        return TEXT("Weak (high µM)");
    else
        return TEXT("Very Weak (mM+)");
}

void AMMGBSA::RunParameterSweep(const TArray<float> &GBScales, const TArray<float> &ChargeScales, const TArray<FString> &LigandKeys)
{
    // Basic parameter sweep: for each GBScale and ChargeScaling, compute binding affinities for the provided ligands (or all ligands if empty)
    if (!ViewerReference)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: RunParameterSweep called but ViewerReference is null."));
        return;
    }

    TArray<FString> Keys;
    if (LigandKeys.Num() > 0)
    {
        Keys = LigandKeys;
    }
    else
    {
        ViewerReference->ClearOverlapMarkers();
        LoadAtomsFromViewer();
        LigandAtoms.GetKeys(Keys);
    }

    if (Keys.Num() == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: No ligand keys provided or found for sweep."));
        return;
    }

    // Save originals
    const float OrigGBScale = GBScale;
    const float OrigChargeScaling = ChargeScaling;

    struct FSweepRow
    {
        float GB;
        float Charge;
        int32 CountValid;
        int32 CountInTarget;
        float MeanDeltaG;
        float MeanAbsDev;
    };
    TArray<FSweepRow> Rows;

    // Target range for 'reasonable' binding energies
    const float TargetLow = -20.0f;
    const float TargetHigh = 0.0f;
    const float TargetMid = (TargetLow + TargetHigh) * 0.5f;

    for (float GB : GBScales)
    {
        for (float CS : ChargeScales)
        {
            GBScale = GB;
            ChargeScaling = CS;

            // Reload atoms so charges are reassigned under new ChargeScaling
            LoadAtomsFromViewer();

            // Also scale receptor atom charges (they were set directly from typical values)
            for (FMMAtom &RA : ReceptorAtoms)
            {
                RA.Charge *= ChargeScaling;
            }

            // Clear cache so each run is fresh
            CachedResults.Empty();

            int32 CountValid = 0;
            int32 CountInTarget = 0;
            double SumDeltaG = 0.0;
            double SumAbsDev = 0.0;

            for (const FString &K : Keys)
            {
                FBindingAffinityResult R = CalculateBindingAffinity(K);
                if (!R.bIsValid)
                    continue;
                CountValid++;
                SumDeltaG += R.DeltaG_Binding;
                double AbsDev = FMath::Abs(R.DeltaG_Binding - TargetMid);
                SumAbsDev += AbsDev;
                if (R.DeltaG_Binding >= TargetLow && R.DeltaG_Binding <= TargetHigh)
                {
                    CountInTarget++;
                }
            }

            float MeanDeltaG = CountValid > 0 ? (float)(SumDeltaG / CountValid) : 0.0f;
            float MeanAbsDev = CountValid > 0 ? (float)(SumAbsDev / CountValid) : FLT_MAX;

            FSweepRow Row;
            Row.GB = GB;
            Row.Charge = CS;
            Row.CountValid = CountValid;
            Row.CountInTarget = CountInTarget;
            Row.MeanDeltaG = MeanDeltaG;
            Row.MeanAbsDev = MeanAbsDev;
            Rows.Add(Row);

            UE_LOG(LogTemp, Log, TEXT("MMGBSA: Sweep GB=%.3f CS=%.3f -> Valid=%d InRange=%d MeanΔG=%.3f MeanAbsDev=%.3f"), GB, CS, CountValid, CountInTarget, MeanDeltaG, MeanAbsDev);
        }
    }

    // Restore originals
    GBScale = OrigGBScale;
    ChargeScaling = OrigChargeScaling;
    CachedResults.Empty();
    LoadAtomsFromViewer();

    // Sort rows: primary by CountInTarget descending, secondary by MeanAbsDev ascending
    Rows.Sort([](const FSweepRow &A, const FSweepRow &B)
              {
        if (A.CountInTarget != B.CountInTarget) return A.CountInTarget > B.CountInTarget;
        return A.MeanAbsDev < B.MeanAbsDev; });

    // Write CSV
    TArray<FString> Lines;
    Lines.Add(TEXT("GBScale,ChargeScaling,CountValid,CountInTarget,MeanDeltaG,MeanAbsDev"));
    for (const FSweepRow &R : Rows)
    {
        Lines.Add(FString::Printf(TEXT("%.6f,%.6f,%d,%d,%.6f,%.6f"), R.GB, R.Charge, R.CountValid, R.CountInTarget, R.MeanDeltaG, R.MeanAbsDev));
    }

    FString FileName = FPaths::Combine(FPaths::ProjectSavedDir(), FString::Printf(TEXT("MMGBSA_Sweep_%s.csv"), *FDateTime::Now().ToString()));
    if (FFileHelper::SaveStringArrayToFile(Lines, *FileName))
    {
        UE_LOG(LogTemp, Log, TEXT("MMGBSA: Parameter sweep saved to %s"), *FileName);
    }
    else
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Failed to write sweep results to %s"), *FileName);
    }
}

bool AMMGBSA::ValidateStructureQuality(const TArray<FMMAtom> &Ligand, FBindingAffinityResult &Result)
{
    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: ValidateStructureQuality - Checking %d ligand atoms against %d receptor atoms"),
           Ligand.Num(), ReceptorAtoms.Num());
    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: MinAllowedDistanceA=%.2f Å, MaxAllowedOverlaps=%d"),
           MinAllowedDistanceA, MaxAllowedOverlaps);

    int32 SevereOverlapCount = 0;
    float MinDistFound = FLT_MAX;

    // Track first few overlaps for detailed logging
    struct FOverlapDetail
    {
        FString ReceptorKey;
        int32 ReceptorIdx;
        FString ReceptorElement;
        int32 LigandIdx;
        FString LigandElement;
        float Distance;
    };
    TArray<FOverlapDetail> FirstOverlaps;

    for (const FMMAtom &R : ReceptorAtoms)
    {
        for (int32 LIdx = 0; LIdx < Ligand.Num(); ++LIdx)
        {
            const FMMAtom &L = Ligand[LIdx];
            float Dist = FVector::Dist(R.Position, L.Position) / 100.0f; // Å
            MinDistFound = FMath::Min(MinDistFound, Dist);

            if (Dist < MinAllowedDistanceA)
            {
                SevereOverlapCount++;

                // Log first 5 overlaps in detail
                if (FirstOverlaps.Num() < 5)
                {
                    FOverlapDetail Detail;
                    Detail.ReceptorKey = R.SourceKey;
                    Detail.ReceptorIdx = R.SourceIndex;
                    Detail.ReceptorElement = R.Element;
                    Detail.LigandIdx = LIdx;
                    Detail.LigandElement = L.Element;
                    Detail.Distance = Dist;
                    FirstOverlaps.Add(Detail);
                }
            }
        }
    }

    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Validation results - Found %d overlaps, MinDist=%.3f Å"),
           SevereOverlapCount, MinDistFound);

    // Log first few overlaps
    for (const FOverlapDetail &Detail : FirstOverlaps)
    {
        UE_LOG(LogTemp, Warning, TEXT("  Overlap: Receptor %s[%d](%s) <-> Ligand[%d](%s) at %.3f Å"),
               *Detail.ReceptorKey, Detail.ReceptorIdx, *Detail.ReceptorElement,
               Detail.LigandIdx, *Detail.LigandElement, Detail.Distance);
    }

    if (SevereOverlapCount > MaxAllowedOverlaps)
    {
        Result.bNeedsMinimization = true;
        Result.QualityWarning = FString::Printf(
            TEXT("Structure has %d severe overlaps (min dist: %.2f Å). Energy minimization required for accurate results."),
            SevereOverlapCount, MinDistFound);

        UE_LOG(LogTemp, Error, TEXT("MMGBSA: VALIDATION FAILED - %s"), *Result.QualityWarning);
        UE_LOG(LogTemp, Error, TEXT("MMGBSA: Threshold: %d overlaps allowed, found %d"), MaxAllowedOverlaps, SevereOverlapCount);

        // For crystal structures with hydrogens, overlaps are common and often acceptable
        // Check if soft-core VDW is enabled - if so, allow calculation
        if (bUseSoftCoreVDW && SevereOverlapCount < 200)
        {
            // Soft-core VDW can handle overlaps by smoothing the potential
            // Only reject if overlaps are catastrophically bad (< 0.5 Å or > 200 overlaps)
            if (MinDistFound < 0.5f)
            {
                UE_LOG(LogTemp, Error, TEXT("MMGBSA: Overlaps are catastrophic (min dist %.2f Å < 0.5 Å). Cannot proceed even with soft-core VDW."), MinDistFound);
                Result.bIsValid = false;
                Result.AffinityClass = TEXT("Invalid - Catastrophic Overlaps");
                return false;
            }

            UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Found %d overlaps (min dist %.2f Å), but soft-core VDW is ENABLED. Proceeding with softened potentials."),
                   SevereOverlapCount, MinDistFound);
            Result.QualityWarning = FString::Printf(TEXT("%d overlaps detected (min %.2f Å), using soft-core VDW"),
                                                    SevereOverlapCount, MinDistFound);
            return true; // Allow calculation to proceed with soft-core
        }

        if (bAutoMinimizeOnOverlap)
        {
            UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Auto-minimization enabled, will attempt to fix structure..."));
            return false;
        }

        Result.bIsValid = false;
        Result.AffinityClass = TEXT("Invalid - Minimize Structure First");
        UE_LOG(LogTemp, Error, TEXT("MMGBSA: Auto-minimize is OFF - calculation will be rejected"));
        return false;
    }

    if (SevereOverlapCount > 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Found %d minor overlaps (min dist: %.2f Å). Results may be less accurate."),
               SevereOverlapCount, MinDistFound);
        Result.QualityWarning = FString::Printf(TEXT("%d minor overlaps detected"), SevereOverlapCount);
    }
    else
    {
        UE_LOG(LogTemp, Log, TEXT("MMGBSA: Structure quality check PASSED - no severe overlaps"));
    }

    return true;
}

bool AMMGBSA::MinimizeStructure(const FString &LigandKey, int32 MaxSteps, float Tolerance)
{
    TArray<FMMAtom> *LigandPtr = LigandAtoms.Find(LigandKey);
    if (!LigandPtr || LigandPtr->Num() == 0)
    {
        UE_LOG(LogTemp, Error, TEXT("MMGBSA: Cannot minimize - ligand %s not found"), *LigandKey);
        return false;
    }

    TArray<FMMAtom> &Ligand = *LigandPtr;

    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Starting minimization for %s (%d steps, tolerance=%.3f)"),
           *LigandKey, MaxSteps, Tolerance);

    MinimizeLigandPosition(Ligand, MaxSteps, Tolerance);

    // Update the viewer's ligand positions if available
    if (ViewerReference)
    {
        FLigandInfo **LigInfoPtr = ViewerReference->LigandMap.Find(LigandKey);
        if (LigInfoPtr && *LigInfoPtr)
        {
            FLigandInfo *LigInfo = *LigInfoPtr;
            for (int32 i = 0; i < Ligand.Num() && i < LigInfo->AtomMeshes.Num(); ++i)
            {
                if (LigInfo->AtomMeshes[i])
                {
                    LigInfo->AtomMeshes[i]->SetWorldLocation(Ligand[i].Position);
                }
            }
        }
    }

    // Clear cached result for this ligand
    CachedResults.Remove(LigandKey);

    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Minimization complete for %s"), *LigandKey);
    return true;
}

void AMMGBSA::MinimizeLigandPosition(TArray<FMMAtom> &Ligand, int32 MaxSteps, float Tolerance)
{
    // Use adaptive step size that decreases as we converge
    float StepSize = 0.5f; // Start with larger steps for severe overlaps
    const float MinStepSize = 0.01f;
    const float StepDecay = 0.98f;

    float PrevEnergy = CalculateMMEnergyPair(ReceptorAtoms, Ligand);
    int32 BadSteps = 0; // Count steps that increase energy

    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Initial energy = %.2f kcal/mol"), PrevEnergy);

    for (int32 Step = 0; Step < MaxSteps; ++Step)
    {
        TArray<FVector> Forces;
        Forces.SetNum(Ligand.Num());

        CalculateForces(ReceptorAtoms, Ligand, Forces);

        // Move atoms along force direction
        float MaxForce = 0.0f;
        for (int32 i = 0; i < Ligand.Num(); ++i)
        {
            float ForceMag = Forces[i].Size();
            MaxForce = FMath::Max(MaxForce, ForceMag);

            if (ForceMag > 0.001f)
            {
                FVector Direction = Forces[i] / ForceMag;
                // Cap individual atom displacement to prevent wild movements
                float Displacement = FMath::Min(StepSize * ForceMag, 1.0f); // Max 1 Å per step
                Ligand[i].Position += Direction * Displacement * 100.0f;    // Convert to UE units
            }
        }

        float CurrentEnergy = CalculateMMEnergyPair(ReceptorAtoms, Ligand);

        // Adaptive step size: reduce if energy increases
        if (CurrentEnergy > PrevEnergy)
        {
            BadSteps++;
            StepSize *= 0.5f; // Halve step size if we're going uphill
            StepSize = FMath::Max(StepSize, MinStepSize);

            if (BadSteps > 10)
            {
                UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Minimization struggling (step %d), energy increased 10 times in a row"), Step);
            }
        }
        else
        {
            BadSteps = 0;
            StepSize *= StepDecay; // Gradually reduce step size as we converge
            StepSize = FMath::Max(StepSize, MinStepSize);
        }

        // Check convergence
        if (MaxForce < Tolerance)
        {
            UE_LOG(LogTemp, Log, TEXT("MMGBSA: Minimization converged at step %d (max force: %.4f, final energy: %.2f)"),
                   Step, MaxForce, CurrentEnergy);
            return;
        }

        // Log progress every 50 steps
        if (Step % 50 == 0)
        {
            UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: Minimize step %d: Energy=%.2f, MaxForce=%.4f, StepSize=%.4f"),
                   Step, CurrentEnergy, MaxForce, StepSize);
        }

        PrevEnergy = CurrentEnergy;
    }

    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Minimization reached max steps (%d) without full convergence"), MaxSteps);
}

void AMMGBSA::CalculateForces(const TArray<FMMAtom> &Receptor, const TArray<FMMAtom> &Ligand, TArray<FVector> &Forces) const
{
    Forces.SetNum(Ligand.Num());
    for (int32 i = 0; i < Forces.Num(); ++i)
    {
        Forces[i] = FVector::ZeroVector;
    }

    const float CutoffA = 8.0f;
    const float MinDist = 1.0f; // Minimum distance to consider (Å)

    for (int32 i = 0; i < Ligand.Num(); ++i)
    {
        const FMMAtom &L = Ligand[i];

        for (const FMMAtom &R : Receptor)
        {
            FVector Diff = L.Position - R.Position;
            float DistUE = Diff.Size();
            float Dist = DistUE / 100.0f; // Å

            // Skip if too far or essentially at same position
            if (Dist < MinDist || Dist > CutoffA)
                continue;

            FVector Direction = Diff / DistUE; // Normalized direction

            // VDW force (repulsive at short range)
            float Epsilon = FMath::Sqrt(GetVDWEpsilon(L.Element) * GetVDWEpsilon(R.Element));
            float Sigma = (GetVDWRadius(L.Element) + GetVDWRadius(R.Element)) * 0.5f;

            // For severe overlaps, use stronger repulsion
            float EffectiveDist = FMath::Max(Dist, 1.2f); // Don't let it go below 1.2 Å
            float SigmaOverR = Sigma / EffectiveDist;
            float SigmaOverR6 = FMath::Pow(SigmaOverR, 6.0f);
            float SigmaOverR7 = SigmaOverR6 * SigmaOverR;
            float SigmaOverR13 = SigmaOverR6 * SigmaOverR7;

            // Lennard-Jones force: F = -dE/dr = 24ε/r * (2(σ/r)^12 - (σ/r)^6)
            float ForceMag = 24.0f * Epsilon * (2.0f * SigmaOverR13 - SigmaOverR7) / EffectiveDist;

            // For very close contacts, add extra repulsion
            if (Dist < 1.5f)
            {
                float ExtraRepulsion = 100.0f * (1.5f - Dist); // Strong push when < 1.5 Å
                ForceMag += ExtraRepulsion;
            }

            // Electrostatic force (usually weaker, less important for overlaps)
            float CoulombK = 332.0636f / 4.0f;
            float ElecForce = CoulombK * L.Charge * R.Charge / (EffectiveDist * EffectiveDist);

            ForceMag += ElecForce;

            // Cap force to prevent explosions
            ForceMag = FMath::Clamp(ForceMag, -200.0f, 200.0f);

            Forces[i] += Direction * ForceMag;
        }
    }

    // Apply damping to total force to prevent oscillations
    for (int32 i = 0; i < Forces.Num(); ++i)
    {
        float ForceMag = Forces[i].Size();
        if (ForceMag > 100.0f)
        {
            // Scale down very large forces
            Forces[i] *= (100.0f / ForceMag);
        }
    }
}

void AMMGBSA::CheckParameterChanges()
{
    bool bNeedsReload = false;

    if (!FMath::IsNearlyEqual(ChargeScaling, LastChargeScaling))
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: ChargeScaling changed from %.3f to %.3f - clearing cache and reloading atoms"),
               LastChargeScaling, ChargeScaling);
        bNeedsReload = true;
        LastChargeScaling = ChargeScaling;
    }

    if (!FMath::IsNearlyEqual(GBScale, LastGBScale))
    {
        UE_LOG(LogTemp, Log, TEXT("MMGBSA: GBScale changed from %.3f to %.3f - clearing cache"),
               LastGBScale, GBScale);
        // GB scaling doesn't affect atom charges, just clear results
        CachedResults.Empty();
        LastGBScale = GBScale;
    }

    if (bNeedsReload)
    {
        CachedResults.Empty();
        LoadAtomsFromViewer();
    }
}
