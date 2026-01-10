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
    TArray<AActor*> Found;
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
    
    auto Add = [this](const FString& Sym, float Rad, float Eps, float Chg, float GBScale)
    {
        FElementData Data;
        Data.VDWRadius = Rad;
        Data.VDWEpsilon = Eps;
        Data.TypicalCharge = Chg;
        Data.GBRadiusScale = GBScale;
        ElementDatabase.Add(Sym, Data);
    };
    
    // Radii in Ångstroms (not scaled), epsilon in kcal/mol
    Add(TEXT("H"),  1.20f, 0.0157f, 0.0f,   0.85f);
    Add(TEXT("C"),  1.70f, 0.0860f, 0.0f,   1.00f);
    Add(TEXT("N"),  1.55f, 0.1700f, -0.3f,  1.05f);
    Add(TEXT("O"),  1.52f, 0.2100f, -0.4f,  1.10f);
    Add(TEXT("S"),  1.80f, 0.2500f, -0.2f,  1.15f);
    Add(TEXT("P"),  1.80f, 0.2000f, 0.5f,   1.10f);
    Add(TEXT("F"),  1.47f, 0.0610f, -0.2f,  1.00f);
    Add(TEXT("CL"), 1.75f, 0.2650f, -0.1f,  1.05f);
    Add(TEXT("BR"), 1.85f, 0.3200f, -0.1f,  1.10f);
    Add(TEXT("I"),  1.98f, 0.4000f, -0.1f,  1.15f);
}

void AMMGBSA::InitializeFromViewer(APDBViewer* Viewer)
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
    if (!ViewerReference) return;
    
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
        for (const auto& Pair : ViewerReference->ResidueMap)
        {
            const FResidueInfo* ResInfo = Pair.Value;
            if (!ResInfo || !ResInfo->bIsVisible) continue;

            for (int32 i = 0; i < ResInfo->AtomPositions.Num(); ++i)
            {
                FMMAtom Atom;
                Atom.Position = ResInfo->AtomPositions[i];
                Atom.bIsReceptor = true;

                if (ResInfo->AtomElements.IsValidIndex(i) && !ResInfo->AtomElements[i].IsEmpty())
                    Atom.Element = ResInfo->AtomElements[i].ToUpper();
                else
                    Atom.Element = TEXT("C");

                const FElementData* Data = ElementDatabase.Find(Atom.Element);
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
                    const FElementData* Def = ElementDatabase.Find(Atom.Element);
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
    for (const auto& Pair : ViewerReference->LigandMap)
    {
        const FString& Key = Pair.Key;
        FLigandInfo* LigInfo = Pair.Value;
        if (!LigInfo) continue;
        
        TArray<FMMAtom> Atoms;
        
            TSet<FString> UnknownElements;
        for (int32 i = 0; i < LigInfo->AtomMeshes.Num(); ++i)
        {
            UStaticMeshComponent* Mesh = LigInfo->AtomMeshes[i];
            if (!Mesh) continue;
            
            FMMAtom Atom;
            Atom.Position = Mesh->GetComponentLocation();
            Atom.bIsReceptor = false;
            
            // Use element information from LigInfo if available
            if (LigInfo->AtomElements.IsValidIndex(i) && !LigInfo->AtomElements[i].IsEmpty())
                Atom.Element = LigInfo->AtomElements[i].ToUpper();
            else
                Atom.Element = TEXT("C"); // Default
            
            const FElementData* Data = ElementDatabase.Find(Atom.Element);
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
                const FElementData* Def = ElementDatabase.Find(Atom.Element);
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
FBindingAffinityResult AMMGBSA::CalculateBindingAffinity(const FString& LigandKey)
{
    FBindingAffinityResult Result;
    Result.LigandKey = LigandKey;
    
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Request to calculate affinity for %s"), *LigandKey);    
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
    TArray<FMMAtom>* LigandPtr = LigandAtoms.Find(LigandKey);
    if (!LigandPtr || LigandPtr->Num() == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Ligand %s not found or has no atoms; ensure the viewer has loaded the molecule and MMGBSA has been initialized"), *LigandKey);
        return Result; // bIsValid remains false
    }
    
    TArray<FMMAtom>& Ligand = *LigandPtr;
    Result.LigandName = LigandKey;
    
    // Create complex (receptor + ligand)
    TArray<FMMAtom> Complex = ReceptorAtoms;
    Complex.Append(Ligand);
    
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Calculating for %s (%d atoms)"), *LigandKey, Ligand.Num());

    // Diagnostic: list ligand atoms with element/position/charge/radius
    float TotalAbsCharge = 0.0f;
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Atoms for %s:"), *LigandKey);
    for (int32 i = 0; i < Ligand.Num(); ++i)
    {
        const FMMAtom& A = Ligand[i];
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
        const FMMAtom& R = ReceptorAtoms[i];
        for (int32 j = 0; j < Ligand.Num(); ++j)
        {
            const FMMAtom& L = Ligand[j];
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
                if (!FMath::IsFinite(PairE)) PairE = FMath::Sign(PairE) * MaxPairVDW;
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

    // Calculate MM energies
    Result.EMM_Complex = CalculateMMEnergy(Complex);
    Result.EMM_Receptor = CalculateMMEnergy(ReceptorAtoms);
    Result.EMM_Ligand = CalculateMMEnergy(Ligand);
    // Use receptor-ligand pairwise energy for ΔEMM to avoid catastrophic cancellation of large internal terms
    Result.DeltaEMM = CalculateMMEnergyPair(ReceptorAtoms, Ligand);
    UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: ΔEMM computed from receptor-ligand pairwise interactions = %.6f"), Result.DeltaEMM);
    
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
    // Use a conservative estimate: ~0.005 kcal/mol per Ų of buried surface
    float BuriedSASA = -(Result.SA_Complex - Result.SA_Receptor - Result.SA_Ligand); // Make positive
    float DesolvationPenalty = BuriedSASA * 0.005f; // Positive penalty (conservative estimate)
    
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
    
    // Convert to Ki
    Result.Ki_uM = DeltaGToKi(Result.DeltaG_Binding);
    Result.AffinityClass = ClassifyAffinity(Result.DeltaG_Binding);
    
    // Log breakdown of components for diagnostics
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: %s | ΔG = %.2f kcal/mol | Ki = %.2f µM | %s"), *LigandKey, Result.DeltaG_Binding, Result.Ki_uM, *Result.AffinityClass);
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Components for %s - ΔEMM=%.3f ΔGB=%.3f ΔSA=%.3f (EMM_C=%.3f EMR=%.3f EML=%.3f)"), *LigandKey, Result.DeltaEMM, Result.DeltaGGB, Result.DeltaGSA, Result.EMM_Complex, Result.EMM_Receptor, Result.EMM_Ligand);

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
    
    for (const auto& Pair : LigandAtoms)
    {
        FBindingAffinityResult Result = CalculateBindingAffinity(Pair.Key);
        Results.Add(Result);
    }
    
    // Sort by binding affinity (most negative = strongest)
    Results.Sort([](const FBindingAffinityResult& A, const FBindingAffinityResult& B)
    {
        return A.DeltaG_Binding < B.DeltaG_Binding;
    });
    
    return Results;
}

bool AMMGBSA::GetCachedResult(const FString& LigandKey, FBindingAffinityResult& OutResult) const
{
    const FBindingAffinityResult* Found = CachedResults.Find(LigandKey);
    if (Found)
    {
        OutResult = *Found;
        return true;
    }
    return false;
}

// MM Energy = Electrostatic + Van der Waals
float AMMGBSA::CalculateMMEnergy(const TArray<FMMAtom>& Atoms) const
{
    float Elec = CalculateElectrostaticEnergy(Atoms);
    float VDW = CalculateVanDerWaalsEnergy(Atoms);
    return Elec + VDW;
}

// Pairwise MM energy between two groups (used for receptor-ligand interactions to avoid large cancellation)
float AMMGBSA::CalculateMMEnergyPair(const TArray<FMMAtom>& Atoms1, const TArray<FMMAtom>& Atoms2) const
{
    float Energy = 0.0f;
    const float CutoffA = this->ReceptorLigandCutoffA; // Å, configurable

    for (int32 i = 0; i < Atoms1.Num(); ++i)
    {
        const FMMAtom& A = Atoms1[i];
        for (int32 j = 0; j < Atoms2.Num(); ++j)
        {
            const FMMAtom& B = Atoms2[j];
            float Dist = FVector::Dist(A.Position, B.Position) / 100.0f; // Å
            if (Dist < 0.5f) Dist = 0.5f;
            if (Dist > CutoffA) continue; // Only consider close receptor-ligand pairs

            // Distance-dependent dielectric to attenuate electrostatics with separation
            float PairDielectric = FMath::Max(4.0f, 4.0f * Dist);
            float CoulombKPair = 332.0636f / PairDielectric;
            float Elec = CoulombKPair * A.Charge * B.Charge / Dist;

            // VDW (Lennard-Jones) with optional soft-core smoothing
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
            if (!FMath::IsFinite(VDW)) VDW = FMath::Sign(VDW) * MaxPairVDW;
            VDW = FMath::Clamp(VDW, -MaxPairVDW, MaxPairVDW);

            // Apply distance-dependent scaling to reduce contribution of very close interactions
            // This helps mitigate the impact of overlapping atoms and makes energies more realistic
            float DistScale = 1.0f;
            if (Dist < 2.0f)
            {
                // Scale down very close interactions (likely overlaps) - less aggressive to allow favorable close contacts
                DistScale = FMath::Lerp(0.5f, 1.0f, FMath::Clamp((Dist - 0.5f) / 1.5f, 0.0f, 1.0f));
            }
            else if (Dist > 4.0f)
            {
                // Scale down distant interactions slightly
                DistScale = FMath::Lerp(1.0f, 0.7f, FMath::Clamp((Dist - 4.0f) / (CutoffA - 4.0f), 0.0f, 1.0f));
            }

            Energy += (Elec + VDW) * DistScale;
        }
    }

    return Energy;
}

float AMMGBSA::CalculateElectrostaticEnergy(const TArray<FMMAtom>& Atoms) const
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
            if (Dist < 0.5f) continue; // Avoid singularity
            if (Dist > CutoffA) continue; // Skip long-range pairs to avoid huge cumulative sums
            
            float E = CoulombK * Atoms[i].Charge * Atoms[j].Charge / Dist;
            Energy += E;
        }
    }
    
    return Energy;
}

float AMMGBSA::CalculateVanDerWaalsEnergy(const TArray<FMMAtom>& Atoms) const
{
    // Lennard-Jones: E = Σᵢⱼ 4εᵢⱼ[(σᵢⱼ/rᵢⱼ)¹² - (σᵢⱼ/rᵢⱼ)⁶]
    float Energy = 0.0f;
    const float CutoffA = 12.0f; // Å - ignore interactions beyond this distance
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        for (int32 j = i + 1; j < Atoms.Num(); ++j)
        {
            float Dist = FVector::Dist(Atoms[i].Position, Atoms[j].Position) / 100.0f; // Convert to Å
            if (Dist < 0.5f) continue;
            if (Dist > CutoffA) continue; // Skip long-range VDW pairs
            
            // Lorentz-Berthelot combining rules
            float Epsilon = FMath::Sqrt(
                GetVDWEpsilon(Atoms[i].Element) * 
                GetVDWEpsilon(Atoms[j].Element)
            );
            
            float Sigma = (GetVDWRadius(Atoms[i].Element) + 
                          GetVDWRadius(Atoms[j].Element)) * 0.5f;

            float DistUse = Dist;
            if (bUseSoftCoreVDW)
            {
                DistUse = FMath::Sqrt(Dist * Dist + SoftCoreAlpha * SoftCoreAlpha);
            }
            
            float SigmaOverR = Sigma / FMath::Max(0.0001f, DistUse);
            float SigmaOverR6 = FMath::Pow(SigmaOverR, 6.0f);
            float SigmaOverR12 = SigmaOverR6 * SigmaOverR6;
            
            float E = 4.0f * Epsilon * (SigmaOverR12 - SigmaOverR6);
            if (!FMath::IsFinite(E)) E = FMath::Sign(E) * MaxPairVDW;
            E = FMath::Clamp(E, -MaxPairVDW, MaxPairVDW);
            Energy += E;
        }
    }
    
    return Energy;
}

// Generalized Born implicit solvation
float AMMGBSA::CalculateGBEnergy(const TArray<FMMAtom>& Atoms)
{
    TArray<FMMAtom> AtomsCopy = Atoms;
    CalculateBornRadii(AtomsCopy);
    const float Factor = -0.5f * 332.0636f * (1.0f / this->InteriorDielectric - 1.0f / this->ExteriorDielectric);
    float Energy = 0.0f;

    // Per-atom self-energy with clamping to avoid huge spikes
    for (const FMMAtom& Atom : AtomsCopy)
    {
        if (Atom.GBRadius > 0.0f)
        {
            float PerAtom = Factor * Atom.Charge * Atom.Charge / Atom.GBRadius;
            if (!FMath::IsFinite(PerAtom)) PerAtom = FMath::Sign(PerAtom) * MaxGBPerAtom;
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
            if (Dist > CutoffA) continue;

            float Ai = AtomsCopy[i].GBRadius;
            float Aj = AtomsCopy[j].GBRadius;
            if (Ai <= 0.0f || Aj <= 0.0f) continue;
            float AiAj = Ai * Aj;
            float Exp = FMath::Exp(-Dist * Dist / (4.0f * AiAj));
            float fGB = 1.0f / FMath::Sqrt(Dist * Dist + AiAj * Exp);

            float PairTerm = 2.0f * Factor * AtomsCopy[i].Charge * AtomsCopy[j].Charge * fGB;
            if (!FMath::IsFinite(PairTerm)) PairTerm = FMath::Sign(PairTerm) * MaxGBPair;
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
float AMMGBSA::CalculateGBEnergyPair(const TArray<FMMAtom>& Atoms1, const TArray<FMMAtom>& Atoms2)
{
    // Create combined system to compute Born radii in binding context (more accurate)
    TArray<FMMAtom> Combined = Atoms1;
    Combined.Append(Atoms2);
    CalculateBornRadii(Combined);
    
    // Extract the updated Born radii back to separate arrays for easier indexing
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
    const float CutoffA = ReceptorLigandGBCutoffA; // Å, configurable
    
    // Distance-dependent scaling to reduce contribution of distant interactions
    // This helps account for the fact that distant interactions are less important for binding
    const float DistScaleStart = 3.0f; // Å - start scaling down beyond this distance
    const float DistScaleEnd = CutoffA; // Å - full scaling at cutoff
    
    // Pairwise GB interactions between the two groups
    for (int32 i = 0; i < Atoms1Copy.Num(); ++i)
    {
        const FMMAtom& A = Atoms1Copy[i];
        for (int32 j = 0; j < Atoms2Copy.Num(); ++j)
        {
            const FMMAtom& B = Atoms2Copy[j];
            float Dist = FVector::Dist(A.Position, B.Position) / 100.0f; // Å
            
            // Only include pairwise GB contribution if within cutoff
            if (Dist > CutoffA || Dist < 0.5f) continue;
            
            float Ai = A.GBRadius;
            float Aj = B.GBRadius;
            if (Ai <= 0.0f || Aj <= 0.0f) continue;
            
            float AiAj = Ai * Aj;
            float Exp = FMath::Exp(-Dist * Dist / (4.0f * AiAj));
            float fGB = 1.0f / FMath::Sqrt(Dist * Dist + AiAj * Exp);
            
            float PairTerm = 2.0f * Factor * A.Charge * B.Charge * fGB;
            
            // Apply distance-dependent scaling to reduce contribution of distant interactions
            float DistScale = 1.0f;
            if (Dist > DistScaleStart)
            {
                // Linearly scale down from 1.0 at DistScaleStart to 0.3 at CutoffA
                float T = (Dist - DistScaleStart) / (DistScaleEnd - DistScaleStart);
                DistScale = FMath::Lerp(1.0f, 0.3f, FMath::Clamp(T, 0.0f, 1.0f));
            }
            PairTerm *= DistScale;
            
            if (!FMath::IsFinite(PairTerm)) PairTerm = FMath::Sign(PairTerm) * MaxGBPair;
            PairTerm = FMath::Clamp(PairTerm, -MaxGBPair, MaxGBPair);
            
            Energy += PairTerm;
        }
    }
    
    return Energy;
}

void AMMGBSA::CalculateBornRadii(TArray<FMMAtom>& Atoms)
{
    // Simplified Born radii calculation
    // More accurate methods: Still, Onufriev-Bashford-Case, etc.
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        float Rho = GetEffectiveBornRadius(Atoms[i], Atoms);
        
        const FElementData* Data = ElementDatabase.Find(Atoms[i].Element);
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

float AMMGBSA::GetEffectiveBornRadius(const FMMAtom& Atom, const TArray<FMMAtom>& AllAtoms)
{
    // Count nearby atoms (burial factor)
    float Rho = 0.0f;
    const float CutoffSq = FMath::Square(1000.0f); // 10 Å
    
    for (const FMMAtom& Other : AllAtoms)
    {
        if (&Other == &Atom) continue;
        
        float DistSq = FVector::DistSquared(Atom.Position, Other.Position);
        if (DistSq < CutoffSq && DistSq > 1.0f)
        {
            float Dist = FMath::Sqrt(DistSq) / 100.0f; // Å
            Rho += FMath::Exp(-Dist / 3.0f); // Exponential decay
        }
    }
    
    return Rho;
}

float AMMGBSA::CalculateSurfaceAreaEnergy(const TArray<FMMAtom>& Atoms)
{
    TArray<FMMAtom> AtomsCopy = Atoms;
    CalculateSASA(AtomsCopy);
    float TotalSASA = 0.0f;
    for (const FMMAtom& Atom : AtomsCopy)
    {
        TotalSASA += Atom.SASA;
    }
    return this->SurfaceTension * TotalSASA;
}

void AMMGBSA::CalculateSASA(TArray<FMMAtom>& Atoms)
{
    // Simplified SASA calculation using sphere overlap
    // More accurate: Lee-Richards algorithm, Shrake-Rupley
    
    const int32 NumSpherePoints = 92; // Fibonacci sphere sampling
    const float ProbeRad = SolventProbeRadius * 100.0f; // Convert to UE units
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        FMMAtom& Atom = Atoms[i];
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
                FMath::Cos(Phi)
            );
            
            // Check if point is accessible (not inside another atom)
            bool bAccessible = true;
            for (int32 j = 0; j < Atoms.Num(); ++j)
            {
                if (i == j) continue;
                
                float Dist = FVector::Dist(TestPoint, Atoms[j].Position);
                if (Dist < (Atoms[j].Radius + ProbeRad))
                {
                    bAccessible = false;
                    break;
                }
            }
            
            if (bAccessible) AccessiblePoints++;
        }
        
        // SASA = 4πr² * (accessible fraction), convert to ų
        float Fraction = (float)AccessiblePoints / NumSpherePoints;
        Atom.SASA = 4.0f * PI * ExtendedRadius * ExtendedRadius * Fraction / 10000.0f; // To ų
    }
}

void AMMGBSA::AssignPartialCharges(TArray<FMMAtom>& Atoms)
{
    // Simple charge assignment based on element
    // Real implementation would use Gasteiger charges or read from file
    TSet<FString> WarnedElements;
    for (FMMAtom& Atom : Atoms)
    {
        Atom.Charge = EstimateChargeFromElement(Atom.Element) * ChargeScaling;

        const FElementData* Data = ElementDatabase.Find(Atom.Element.ToUpper());
        if (FMath::IsNearlyZero(Atom.Charge, 1e-6f) && !WarnedElements.Contains(Atom.Element))
        {
            if (Data)
            {
                // Common elements (e.g., carbon) may intentionally have 0.0 typical charge; log at Verbose instead of Warning
                if (FMath::IsNearlyZero(Data->TypicalCharge * ChargeScaling, 1e-6f))
                {
                    UE_LOG(LogTemp, Verbose, TEXT("MMGBSA: Assigned zero partial charge for element '%s' (common default after scaling). Consider using a dedicated charge model for more accurate electrostatics."), *Atom.Element);
                }
                else
                {
                    UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Assigned zero partial charge for element '%s' though element has a non-zero typical charge after scaling."), *Atom.Element);
                }
            }
            else
            {
                UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Assigned zero partial charge for unknown element '%s' (defaulting to 0.0)."), *Atom.Element);
            }

            WarnedElements.Add(Atom.Element);
        }
    }
}

float AMMGBSA::EstimateChargeFromElement(const FString& Element) const
{
    const FElementData* Data = ElementDatabase.Find(Element.ToUpper());
    return Data ? Data->TypicalCharge : 0.0f;
}

float AMMGBSA::GetVDWRadius(const FString& Element) const
{
    const FElementData* Data = ElementDatabase.Find(Element.ToUpper());
    return Data ? Data->VDWRadius : 1.70f;
}

float AMMGBSA::GetVDWEpsilon(const FString& Element) const
{
    const FElementData* Data = ElementDatabase.Find(Element.ToUpper());
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
        if (!FMath::IsFinite((float)Ki_M)) bOverflow = true;
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

void AMMGBSA::RunParameterSweep(const TArray<float>& GBScales, const TArray<float>& ChargeScales, const TArray<FString>& LigandKeys)
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

    struct FSweepRow { float GB; float Charge; int32 CountValid; int32 CountInTarget; float MeanDeltaG; float MeanAbsDev; };
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
            for (FMMAtom& RA : ReceptorAtoms)
            {
                RA.Charge *= ChargeScaling;
            }

            // Clear cache so each run is fresh
            CachedResults.Empty();

            int32 CountValid = 0;
            int32 CountInTarget = 0;
            double SumDeltaG = 0.0;
            double SumAbsDev = 0.0;

            for (const FString& K : Keys)
            {
                FBindingAffinityResult R = CalculateBindingAffinity(K);
                if (!R.bIsValid) continue;
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
    Rows.Sort([](const FSweepRow& A, const FSweepRow& B)
    {
        if (A.CountInTarget != B.CountInTarget) return A.CountInTarget > B.CountInTarget;
        return A.MeanAbsDev < B.MeanAbsDev;
    });

    // Write CSV
    TArray<FString> Lines;
    Lines.Add(TEXT("GBScale,ChargeScaling,CountValid,CountInTarget,MeanDeltaG,MeanAbsDev"));
    for (const FSweepRow& R : Rows)
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

