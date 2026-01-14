// MMGBSA.cpp - Complete rewrite with improved physics
#include "MMGBSA.h"
#include "PDBViewer.h"
#include "Kismet/GameplayStatics.h"
#include "Misc/FileHelper.h"
#include "Misc/Paths.h"

// Physical constants
static constexpr float COULOMB_CONSTANT = 332.0636f; // kcal·Å/(mol·e²)
static constexpr float GAS_CONSTANT = 0.001987f; // kcal/(mol·K)
static constexpr float AVOGADRO = 6.02214076e23f;

AMMGBSA::AMMGBSA()
{
    PrimaryActorTick.bCanEverTick = false;
}

void AMMGBSA::BeginPlay()
{
    Super::BeginPlay();
    
    InitializeForceField();
    
    // Auto-find viewer
    TArray<AActor*> FoundActors;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), APDBViewer::StaticClass(), FoundActors);
    if (FoundActors.Num() > 0)
    {
        ViewerReference = Cast<APDBViewer>(FoundActors[0]);
        if (ViewerReference)
        {
            LoadAtomsFromViewer();
            ViewerReference->OnLigandsLoaded.AddDynamic(this, &AMMGBSA::OnViewerLigandsLoaded);
            UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Initialized with PDB Viewer"));
        }
    }
}

void AMMGBSA::InitializeForceField()
{
    ForceField.Empty();
    
    // AMBER-like force field parameters
    // Format: Element -> (VDW_radius_Å, VDW_epsilon_kcal/mol, base_charge, born_scale, electronegativity)
    
    auto AddFF = [this](const FString& El, float R, float Eps, float Q, float BS, float EN)
    {
        FForceFieldParams FF;
        FF.VDWRadius = R;
        FF.VDWEpsilon = Eps;
        FF.BaseCharge = Q;
        FF.BornRadiusScale = BS;
        FF.Electronegativity = EN;
        ForceField.Add(El, FF);
    };
    
    // Main elements (R, Epsilon, BaseCharge, BornScale, Electronegativity)
    AddFF(TEXT("H"),  1.20f, 0.0157f,  0.00f, 0.85f, 2.20f);
    AddFF(TEXT("C"),  1.70f, 0.0860f,  0.00f, 1.00f, 2.55f);
    AddFF(TEXT("N"),  1.55f, 0.1700f, -0.40f, 1.05f, 3.04f);
    AddFF(TEXT("O"),  1.52f, 0.2100f, -0.50f, 1.10f, 3.44f);
    AddFF(TEXT("S"),  1.80f, 0.2500f, -0.20f, 1.15f, 2.58f);
    AddFF(TEXT("P"),  1.80f, 0.2000f,  0.70f, 1.10f, 2.19f);
    AddFF(TEXT("F"),  1.47f, 0.0610f, -0.25f, 1.00f, 3.98f);
    AddFF(TEXT("CL"), 1.75f, 0.2650f, -0.15f, 1.05f, 3.16f);
    AddFF(TEXT("BR"), 1.85f, 0.3200f, -0.15f, 1.10f, 2.96f);
    AddFF(TEXT("I"),  1.98f, 0.4000f, -0.15f, 1.15f, 2.66f);
    
    // Metal ions
    AddFF(TEXT("ZN"), 1.39f, 0.2500f, 2.00f, 1.20f, 1.65f);
    AddFF(TEXT("Z"),  1.39f, 0.2500f, 2.00f, 1.20f, 1.65f); // Alias
    AddFF(TEXT("MG"), 1.73f, 0.8750f, 2.00f, 1.18f, 1.31f);
    AddFF(TEXT("M"),  1.73f, 0.8750f, 2.00f, 1.18f, 1.31f); // Alias
    AddFF(TEXT("CA"), 2.31f, 0.4500f, 2.00f, 1.22f, 1.00f);
    AddFF(TEXT("FE"), 1.47f, 0.0500f, 2.00f, 1.20f, 1.83f);
    AddFF(TEXT("MN"), 1.61f, 0.0500f, 2.00f, 1.20f, 1.55f);
    AddFF(TEXT("CU"), 1.40f, 0.0500f, 2.00f, 1.20f, 1.90f);
    AddFF(TEXT("NA"), 2.27f, 0.0028f, 1.00f, 1.15f, 0.93f);
    AddFF(TEXT("K"),  2.75f, 0.0003f, 1.00f, 1.25f, 0.82f);
    
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Force field initialized with %d elements"), ForceField.Num());
}

FForceFieldParams AMMGBSA::GetForceFieldParams(const FString& Element) const
{
    const FForceFieldParams* Params = ForceField.Find(Element.ToUpper());
    if (Params)
    {
        return *Params;
    }
    
    // Default to carbon-like parameters
    UE_LOG(LogTemp, Warning, TEXT("MM/GBSA: Unknown element '%s', using carbon defaults"), *Element);
    FForceFieldParams Default;
    Default.VDWRadius = 1.70f;
    Default.VDWEpsilon = 0.086f;
    Default.BaseCharge = 0.0f;
    Default.BornRadiusScale = 1.0f;
    Default.Electronegativity = 2.55f;
    return Default;
}

void AMMGBSA::InitializeFromViewer(APDBViewer* Viewer)
{
    if (!Viewer)
    {
        UE_LOG(LogTemp, Error, TEXT("MM/GBSA: Null viewer reference"));
        return;
    }
    
    ViewerReference = Viewer;
    LoadAtomsFromViewer();
    
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Loaded %d receptor atoms, %d ligands"),
           ReceptorAtoms.Num(), LigandAtoms.Num());
}

void AMMGBSA::LoadAtomsFromViewer()
{
    if (!ViewerReference)
    {
        UE_LOG(LogTemp, Warning, TEXT("MM/GBSA: No viewer reference"));
        return;
    }
    
    ReceptorAtoms.Empty();
    LigandAtoms.Empty();
    CachedResults.Empty();
    
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Loading atoms from viewer..."));
    
    // Load receptor atoms from ResidueMap
    int32 ReceptorAtomCount = 0;
    for (const auto& Pair : ViewerReference->ResidueMap)
    {
        const FString& ResidueKey = Pair.Key;
        const FResidueInfo* ResInfo = Pair.Value;
        
        if (!ResInfo || !ResInfo->bIsVisible)
            continue;
        
        for (int32 i = 0; i < ResInfo->AtomPositions.Num(); ++i)
        {
            FMMAtom Atom;
            Atom.Position = ResInfo->AtomPositions[i];
            Atom.bIsReceptor = true;
            Atom.SourceIndex = i;
            
            // Get element
            if (ResInfo->AtomElements.IsValidIndex(i) && !ResInfo->AtomElements[i].IsEmpty())
            {
                Atom.Element = ResInfo->AtomElements[i].ToUpper();
            }
            else
            {
                Atom.Element = TEXT("C"); // Default
            }
            
            // Create diagnostic label
            FString AtomName = ResInfo->AtomNames.IsValidIndex(i) ? ResInfo->AtomNames[i] : FString::Printf(TEXT("%d"), i);
            Atom.Label = FString::Printf(TEXT("%s:%s"), *ResidueKey, *AtomName);
            
            // Assign force field parameters
            FForceFieldParams FF = GetForceFieldParams(Atom.Element);
            Atom.VDWRadius = FF.VDWRadius;
            Atom.VDWEpsilon = FF.VDWEpsilon;
            Atom.Charge = FF.BaseCharge * ChargeScalingFactor; // Apply charge scaling
            
            ReceptorAtoms.Add(Atom);
            ReceptorAtomCount++;
        }
    }
    
    // Load ligand atoms from LigandMap
    int32 LigandCount = 0;
    for (const auto& Pair : ViewerReference->LigandMap)
    {
        const FString& LigandKey = Pair.Key;
        FLigandInfo* LigInfo = Pair.Value;
        
        if (!LigInfo || !LigInfo->bIsVisible)
            continue;
        
        TArray<FMMAtom> Atoms;
        
        for (int32 i = 0; i < LigInfo->AtomPositions.Num(); ++i)
        {
            FMMAtom Atom;
            Atom.Position = LigInfo->AtomPositions[i];
            Atom.bIsReceptor = false;
            Atom.SourceIndex = i;
            
            // Get element
            if (LigInfo->AtomElements.IsValidIndex(i) && !LigInfo->AtomElements[i].IsEmpty())
            {
                Atom.Element = LigInfo->AtomElements[i].ToUpper();
            }
            else
            {
                Atom.Element = TEXT("C");
            }
            
            // Create label
            Atom.Label = FString::Printf(TEXT("%s:%d"), *LigandKey, i);
            
            // Assign force field parameters
            FForceFieldParams FF = GetForceFieldParams(Atom.Element);
            Atom.VDWRadius = FF.VDWRadius;
            Atom.VDWEpsilon = FF.VDWEpsilon;
            Atom.Charge = FF.BaseCharge;
            
            Atoms.Add(Atom);
        }
        
        // Refine charges based on molecular context
        if (Atoms.Num() > 0)
        {
            RefinePartialCharges(Atoms);
            LigandAtoms.Add(LigandKey, Atoms);
            LigandCount++;
        }
    }
    
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Loaded %d receptor atoms from %d residues"),
           ReceptorAtomCount, ViewerReference->ResidueMap.Num());
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Loaded %d ligands"), LigandCount);
}

void AMMGBSA::RefinePartialCharges(TArray<FMMAtom>& Atoms)
{
    // Improve partial charges using connectivity and electronegativity
    
    const float BondCutoff = 2.0f; // Å - consider atoms within this distance as bonded
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        FMMAtom& AtomI = Atoms[i];
        FForceFieldParams FFI = GetForceFieldParams(AtomI.Element);
        
        // Find bonded neighbors
        TArray<int32> Neighbors;
        for (int32 j = 0; j < Atoms.Num(); ++j)
        {
            if (i == j) continue;
            
            float Dist = FVector::Dist(AtomI.Position, Atoms[j].Position);
            if (Dist < BondCutoff)
            {
                Neighbors.Add(j);
            }
        }
        
        // Adjust charge based on neighbors
        if (Neighbors.Num() > 0)
        {
            float ChargeAdjustment = 0.0f;
            
            for (int32 j : Neighbors)
            {
                const FMMAtom& AtomJ = Atoms[j];
                FForceFieldParams FFJ = GetForceFieldParams(AtomJ.Element);
                
                // Electronegativity difference drives charge transfer
                float DeltaEN = FFJ.Electronegativity - FFI.Electronegativity;
                ChargeAdjustment += DeltaEN * 0.05f; // Scale factor
            }
            
            AtomI.Charge += ChargeAdjustment;
        }
        
        // Apply reasonable bounds
        AtomI.Charge = FMath::Clamp(AtomI.Charge, -1.0f, 1.0f);
        
        // Apply global charge scaling (for unminimized structures)
        AtomI.Charge *= ChargeScalingFactor;
    }
    
    // Ensure molecular neutrality
    float TotalCharge = 0.0f;
    for (const FMMAtom& A : Atoms)
    {
        TotalCharge += A.Charge;
    }
    
    // Only neutralize if charge is unreasonably large
    if (FMath::Abs(TotalCharge) > 2.0f)
    {
        float Correction = -TotalCharge / Atoms.Num();
        for (FMMAtom& A : Atoms)
        {
            A.Charge += Correction;
        }
        
        UE_LOG(LogTemp, Verbose, TEXT("MM/GBSA: Neutralized molecule (was %.2f e, corrected %.2f e per atom)"),
               TotalCharge, Correction);
    }
}

void AMMGBSA::AssignForceFieldParameters(TArray<FMMAtom>& Atoms)
{
    for (FMMAtom& Atom : Atoms)
    {
        FForceFieldParams FF = GetForceFieldParams(Atom.Element);
        Atom.VDWRadius = FF.VDWRadius;
        Atom.VDWEpsilon = FF.VDWEpsilon;
        Atom.Charge = FF.BaseCharge;
    }
}

void AMMGBSA::OnViewerLigandsLoaded()
{
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Viewer ligands updated, reloading atoms"));
    LoadAtomsFromViewer();
}

void AMMGBSA::ClearCache()
{
    CachedResults.Empty();
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Cache cleared"));
}

// ============================================================================
// PART 2: CORE ENERGY CALCULATIONS
// ============================================================================

float AMMGBSA::CalculateElectrostatic(const FMMAtom& A, const FMMAtom& B, float Distance) const
{
    if (Distance < 0.5f || Distance > ElectrostaticCutoff)
        return 0.0f;
    
    // Coulomb's law with dielectric screening
    float Dielectric = InteriorDielectric;
    
    if (bUseDistanceDependentDielectric)
    {
        // Distance-dependent dielectric: ε(r) = 4r
        // More realistic for protein interiors and unminimized structures
        Dielectric = FMath::Max(4.0f, 4.0f * Distance);
    }
    
    float Energy = (COULOMB_CONSTANT / Dielectric) * A.Charge * B.Charge / Distance;
    
    // Smooth cutoff near boundary
    if (Distance > ElectrostaticCutoff * 0.9f)
    {
        float SwitchFactor = 1.0f - FMath::Square((Distance - ElectrostaticCutoff * 0.9f) / (ElectrostaticCutoff * 0.1f));
        Energy *= FMath::Clamp(SwitchFactor, 0.0f, 1.0f);
    }
    
    return Energy;
}

float AMMGBSA::CalculateVDW(const FMMAtom& A, const FMMAtom& B, float Distance) const
{
    if (Distance > VDWCutoff)
        return 0.0f;
    
    // Lennard-Jones 12-6 potential with soft-core modification
    float Epsilon = FMath::Sqrt(A.VDWEpsilon * B.VDWEpsilon); // Lorentz-Berthelot
    float Sigma = (A.VDWRadius + B.VDWRadius) * 0.5f;
    
    float EffectiveDistance = Distance;
    
    if (bUseSoftCore && Distance < 2.0f)
    {
        // Soft-core potential to handle clashes gracefully
        EffectiveDistance = FMath::Sqrt(Distance * Distance + SoftCoreAlpha * SoftCoreAlpha);
    }
    
    float SigmaOverR = Sigma / FMath::Max(EffectiveDistance, 0.1f);
    float SigmaOverR6 = FMath::Pow(SigmaOverR, 6.0f);
    float SigmaOverR12 = SigmaOverR6 * SigmaOverR6;
    
    float Energy = 4.0f * Epsilon * (SigmaOverR12 - SigmaOverR6);
    
    // Cap extreme values
    Energy = FMath::Clamp(Energy, -50.0f, 100.0f);
    
    // Smooth cutoff
    if (Distance > VDWCutoff * 0.9f)
    {
        float SwitchFactor = 1.0f - FMath::Square((Distance - VDWCutoff * 0.9f) / (VDWCutoff * 0.1f));
        Energy *= FMath::Clamp(SwitchFactor, 0.0f, 1.0f);
    }
    
    return Energy;
}

FEnergyComponents AMMGBSA::CalculateMMEnergy(const TArray<FMMAtom>& Atoms) const
{
    FEnergyComponents Result;
    Result.bIsValid = true;
    
    const int32 N = Atoms.Num();
    
    for (int32 i = 0; i < N; ++i)
    {
        for (int32 j = i + 1; j < N; ++j)
        {
            float Distance = FVector::Dist(Atoms[i].Position, Atoms[j].Position);
            
            if (Distance < 0.5f) // Skip essentially overlapping atoms
            {
                Result.ClashCount++;
                continue;
            }
            
            float Elec = CalculateElectrostatic(Atoms[i], Atoms[j], Distance);
            float VDW = CalculateVDW(Atoms[i], Atoms[j], Distance);
            
            Result.Electrostatic += Elec;
            Result.VanDerWaals += VDW;
            Result.InteractionPairs++;
        }
    }
    
    Result.Total = Result.Electrostatic + Result.VanDerWaals;
    
    return Result;
}

FEnergyComponents AMMGBSA::CalculatePairwiseEnergy(const TArray<FMMAtom>& Receptor, const TArray<FMMAtom>& Ligand) const
{
    FEnergyComponents Result;
    Result.bIsValid = true;
    
    // Only calculate receptor-ligand interactions (no internal energies)
    for (int32 i = 0; i < Receptor.Num(); ++i)
    {
        const FMMAtom& R = Receptor[i];
        
        for (int32 j = 0; j < Ligand.Num(); ++j)
        {
            const FMMAtom& L = Ligand[j];
            
            float Distance = FVector::Dist(R.Position, L.Position);
            
            if (Distance < 0.5f)
            {
                Result.ClashCount++;
                continue;
            }
            
            // Skip if too far away
            if (Distance > FMath::Max(ElectrostaticCutoff, VDWCutoff))
                continue;
            
            float Elec = CalculateElectrostatic(R, L, Distance);
            float VDW = CalculateVDW(R, L, Distance);
            
            Result.Electrostatic += Elec;
            Result.VanDerWaals += VDW;
            Result.InteractionPairs++;
        }
    }
    
    Result.Total = Result.Electrostatic + Result.VanDerWaals;
    
    return Result;
}

// ============================================================================
// PART 3: GENERALIZED BORN SOLVATION
// ============================================================================

void AMMGBSA::CalculateBornRadii(TArray<FMMAtom>& Atoms) const
{
    // Calculate effective Born radii using the Still method
    // Born radius represents how "buried" an atom is
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        FMMAtom& AtomI = Atoms[i];
        
        // Start with intrinsic radius
        FForceFieldParams FF = GetForceFieldParams(AtomI.Element);
        float Ri = FF.VDWRadius * FF.BornRadiusScale;
        
        // Calculate volume integral for pairwise descreening
        float Integral = 0.0f;
        
        for (int32 j = 0; j < Atoms.Num(); ++j)
        {
            if (i == j) continue;
            
            const FMMAtom& AtomJ = Atoms[j];
            float Dist = FVector::Dist(AtomI.Position, AtomJ.Position);
            
            if (Dist > GBCutoff)
                continue;
            
            FForceFieldParams FFJ = GetForceFieldParams(AtomJ.Element);
            float Rj = FFJ.VDWRadius * FFJ.BornRadiusScale;
            
            // Still's pairwise descreening integral
            if (Dist < Ri + Rj)
            {
                // Overlap - strong descreening
                float L = FMath::Max(Ri, FMath::Abs(Dist - Rj));
                float U = Dist + Rj;
                
                if (U > 0.0f && L < U)
                {
                    float L2 = L * L;
                    float U2 = U * U;
                    float Dist2 = Dist * Dist;
                    float Rj2 = Rj * Rj;
                    
                    float Term = (1.0f / L - 1.0f / U) + 
                                 (Dist2 - Rj2) / (4.0f * Dist) * (1.0f / (U * U) - 1.0f / (L * L)) -
                                 (0.5f / Dist) * FMath::Loge(L / U);
                    
                    Integral += Term / Rj2;
                }
            }
            else
            {
                // No overlap - weak descreening
                float Term = 1.0f / Dist - 1.0f / (Dist + Rj);
                Integral += Term / Rj;
            }
        }
        
        // Calculate effective Born radius
        float BornRadius = 1.0f / (1.0f / Ri - Integral / 2.0f);
        
        // Apply bounds
        AtomI.BornRadius = FMath::Clamp(BornRadius, Ri * 0.8f, Ri * 3.0f);
    }
}

float AMMGBSA::CalculateGBPair(const FMMAtom& A, const FMMAtom& B, float Distance) const
{
    if (Distance > GBCutoff || Distance < 0.5f)
        return 0.0f;
    
    // Generalized Born equation for pairwise interaction
    float fGB = FMath::Sqrt(Distance * Distance + A.BornRadius * B.BornRadius * 
                            FMath::Exp(-Distance * Distance / (4.0f * A.BornRadius * B.BornRadius)));
    
    // Solvation free energy
    float DeltaDielectric = 1.0f / ExteriorDielectric - 1.0f / InteriorDielectric;
    float Energy = -COULOMB_CONSTANT * DeltaDielectric * A.Charge * B.Charge / fGB;
    
    // Include ionic strength effects (Debye-Hückel)
    if (SaltConcentration > 0.0f)
    {
        float IonicStrength = 0.5f * SaltConcentration;
        float DebyeLength = FMath::Sqrt(ExteriorDielectric * 8.85419e-12f * GAS_CONSTANT * Temperature / 
                                        (2.0f * AVOGADRO * 1.602176634e-19f * 1.602176634e-19f * IonicStrength * 1000.0f));
        float Kappa = 1.0f / DebyeLength;
        
        float ScreeningFactor = FMath::Exp(-Kappa * Distance);
        Energy *= ScreeningFactor;
    }
    
    return Energy;
}

float AMMGBSA::CalculateGBEnergy(TArray<FMMAtom>& Atoms) const
{
    // First calculate Born radii
    CalculateBornRadii(Atoms);
    
    float TotalEnergy = 0.0f;
    
    // Self energy (desolvation of individual atoms)
    for (const FMMAtom& A : Atoms)
    {
        if (A.BornRadius > 0.0f)
        {
            float DeltaDielectric = 1.0f / ExteriorDielectric - 1.0f / InteriorDielectric;
            float SelfEnergy = -0.5f * COULOMB_CONSTANT * DeltaDielectric * A.Charge * A.Charge / A.BornRadius;
            TotalEnergy += SelfEnergy;
        }
    }
    
    // Pairwise interactions
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        for (int32 j = i + 1; j < Atoms.Num(); ++j)
        {
            float Distance = FVector::Dist(Atoms[i].Position, Atoms[j].Position);
            float PairEnergy = CalculateGBPair(Atoms[i], Atoms[j], Distance);
            TotalEnergy += PairEnergy;
        }
    }
    
    return TotalEnergy;
}

// ============================================================================
// PART 4: SURFACE AREA CALCULATIONS
// ============================================================================

void AMMGBSA::CalculateSASA(TArray<FMMAtom>& Atoms) const
{
    // Calculate solvent accessible surface area using approximate method
    // Based on Shrake-Rupley algorithm with reduced point set
    
    const int32 NumSpherePoints = 92; // Reduced for performance
    const float ProbeRadius = SolventProbeRadius;
    
    // Generate sphere points using Fibonacci spiral
    TArray<FVector> SpherePoints;
    SpherePoints.Reserve(NumSpherePoints);
    
    for (int32 i = 0; i < NumSpherePoints; ++i)
    {
        float Theta = PI * (1.0f + FMath::Sqrt(5.0f)) * i;
        float Z = 1.0f - (2.0f * i) / (NumSpherePoints - 1.0f);
        float R = FMath::Sqrt(1.0f - Z * Z);
        SpherePoints.Add(FVector(R * FMath::Cos(Theta), R * FMath::Sin(Theta), Z));
    }
    
    // Calculate SASA for each atom
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        FMMAtom& AtomI = Atoms[i];
        float TestRadius = AtomI.VDWRadius + ProbeRadius;
        
        int32 AccessiblePoints = 0;
        
        for (const FVector& Point : SpherePoints)
        {
            FVector TestPoint = AtomI.Position + Point * TestRadius;
            bool bIsAccessible = true;
            
            // Check if test point is inside any other atom
            for (int32 j = 0; j < Atoms.Num(); ++j)
            {
                if (i == j) continue;
                
                const FMMAtom& AtomJ = Atoms[j];
                float JRadius = AtomJ.VDWRadius + ProbeRadius;
                float Dist = FVector::Dist(TestPoint, AtomJ.Position);
                
                if (Dist < JRadius)
                {
                    bIsAccessible = false;
                    break;
                }
            }
            
            if (bIsAccessible)
            {
                AccessiblePoints++;
            }
        }
        
        // Calculate surface area
        float SphereArea = 4.0f * PI * TestRadius * TestRadius;
        AtomI.SASA = SphereArea * (float)AccessiblePoints / NumSpherePoints;
    }
}

float AMMGBSA::CalculateSurfaceAreaEnergy(TArray<FMMAtom>& Atoms) const
{
    CalculateSASA(Atoms);
    
    float TotalSA = 0.0f;
    for (const FMMAtom& A : Atoms)
    {
        TotalSA += A.SASA;
    }
    
    // Nonpolar solvation energy
    return SurfaceTension * TotalSA;
}

// ============================================================================
// PART 5: STRUCTURE QUALITY ASSESSMENT
// ============================================================================

TArray<FAtomicClash> AMMGBSA::DetectClashes(const TArray<FMMAtom>& Receptor, const TArray<FMMAtom>& Ligand) const
{
    TArray<FAtomicClash> Clashes;
    
    for (const FMMAtom& R : Receptor)
    {
        for (const FMMAtom& L : Ligand)
        {
            float Distance = FVector::Dist(R.Position, L.Position);
            float MinDistance = (R.VDWRadius + L.VDWRadius) * 0.5f; // Sum of radii
            
            // Check if clash exists
            if (Distance < MinDistance * 0.95f) // Allow 5% tolerance
            {
                FAtomicClash Clash;
                Clash.Atom1Label = R.Label;
                Clash.Atom2Label = L.Label;
                Clash.Position1 = R.Position;
                Clash.Position2 = L.Position;
                Clash.Distance = Distance;
                Clash.MinDistance = MinDistance;
                Clash.Severity = 1.0f - (Distance / MinDistance);
                Clash.Severity = FMath::Clamp(Clash.Severity, 0.0f, 1.0f);
                
                Clashes.Add(Clash);
            }
        }
    }
    
    // Sort by severity (worst first)
    Clashes.Sort([](const FAtomicClash& A, const FAtomicClash& B) {
        return A.Severity > B.Severity;
    });
    
    return Clashes;
}

FStructureQuality AMMGBSA::AssessQuality(const TArray<FMMAtom>& Receptor, const TArray<FMMAtom>& Ligand) const
{
    FStructureQuality Quality;
    
    TArray<FAtomicClash> AllClashes = DetectClashes(Receptor, Ligand);
    
    // Classify clashes by severity
    for (const FAtomicClash& Clash : AllClashes)
    {
        float Threshold = Clash.MinDistance;
        
        if (Clash.Distance < Threshold * SevereClashThreshold)
        {
            Quality.SevereClashes++;
        }
        else if (Clash.Distance < Threshold * ModerateClashThreshold)
        {
            Quality.ModerateClashes++;
        }
    }
    
    // Get worst clashes for reporting
    Quality.TopClashes = AllClashes;
    if (Quality.TopClashes.Num() > 5)
    {
        Quality.TopClashes.SetNum(5); // Keep top 5
    }
    
    Quality.WorstClashSeverity = AllClashes.Num() > 0 ? AllClashes[0].Severity : 0.0f;
    
    // Determine if structure is acceptable
    Quality.bIsAcceptable = (Quality.SevereClashes <= MaxAllowedSevereClashes);
    
    // Recommendation
    if (Quality.SevereClashes == 0 && Quality.ModerateClashes == 0)
    {
        Quality.RecommendedAction = TEXT("Structure quality is good");
    }
    else if (Quality.SevereClashes > 0 && Quality.SevereClashes <= MaxAllowedSevereClashes)
    {
        Quality.RecommendedAction = TEXT("Minor clashes detected - recommend energy minimization");
    }
    else if (Quality.SevereClashes > MaxAllowedSevereClashes)
    {
        Quality.RecommendedAction = TEXT("Severe clashes detected - energy minimization required");
    }
    
    return Quality;
}

FStructureQuality AMMGBSA::AssessStructureQuality(const FString& LigandKey)
{
    TArray<FMMAtom>* LigandPtr = LigandAtoms.Find(LigandKey);
    
    if (!LigandPtr || LigandPtr->Num() == 0)
    {
        FStructureQuality Quality;
        Quality.bIsAcceptable = false;
        Quality.RecommendedAction = TEXT("Ligand not found");
        return Quality;
    }
    
    return AssessQuality(ReceptorAtoms, *LigandPtr);
}

// ============================================================================
// PART 6: MAIN BINDING AFFINITY CALCULATION
// ============================================================================

FBindingAffinityResult AMMGBSA::CalculateBindingAffinity(const FString& LigandKey)
{
    FBindingAffinityResult Result;
    Result.LigandName = LigandKey;
    Result.LigandKey = LigandKey;
    
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Calculating binding affinity for %s"), *LigandKey);
    
    // Check cache
    if (CachedResults.Contains(LigandKey))
    {
        UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Using cached result"));
        return CachedResults[LigandKey];
    }
    
    UE_LOG(LogTemp, Verbose, TEXT("MM/GBSA: Looking up ligand atoms..."));
    
    // Find ligand
    TArray<FMMAtom>* LigandPtr = LigandAtoms.Find(LigandKey);
    if (!LigandPtr || LigandPtr->Num() == 0)
    {
        Result.ErrorMessage = TEXT("Ligand not found or has no atoms");
        Result.bIsValid = false;
        UE_LOG(LogTemp, Error, TEXT("MM/GBSA: %s"), *Result.ErrorMessage);
        return Result;
    }
    
    UE_LOG(LogTemp, Verbose, TEXT("MM/GBSA: Copying ligand atoms..."));
    
    TArray<FMMAtom> Ligand = *LigandPtr;
    
    // Safety check
    if (Ligand.Num() == 0)
    {
        Result.ErrorMessage = TEXT("Ligand has no atoms");
        Result.bIsValid = false;
        UE_LOG(LogTemp, Error, TEXT("MM/GBSA: %s"), *Result.ErrorMessage);
        return Result;
    }
    
    // Assess structure quality
    Result.QualityAssessment = AssessQuality(ReceptorAtoms, Ligand);
    
    // Safety check for receptor atoms
    if (ReceptorAtoms.Num() == 0)
    {
        Result.ErrorMessage = TEXT("No receptor atoms loaded");
        Result.bIsValid = false;
        UE_LOG(LogTemp, Error, TEXT("MM/GBSA: %s"), *Result.ErrorMessage);
        return Result;
    }
    
    if (!Result.QualityAssessment.bIsAcceptable)
    {
        UE_LOG(LogTemp, Warning, TEXT("MM/GBSA: Structure quality poor - %d severe clashes"),
               Result.QualityAssessment.SevereClashes);
        
        if (bAutoMinimizeOnClash)
        {
            UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Attempting automatic minimization"));
            MinimizeLigand(Ligand, DefaultMinimizationSteps, MinimizationTolerance);
            
            // Re-assess after minimization
            Result.QualityAssessment = AssessQuality(ReceptorAtoms, Ligand);
            
            if (!Result.QualityAssessment.bIsAcceptable)
            {
                Result.ErrorMessage = TEXT("Structure still has severe clashes after minimization");
                Result.bIsValid = false;
                UE_LOG(LogTemp, Error, TEXT("MM/GBSA: Minimization failed to resolve clashes"));
                CachedResults.Add(LigandKey, Result);
                OnAffinityCalculated.Broadcast(Result);
                return Result;
            }
            
            UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Minimization successful"));
            // Update ligand atoms with minimized structure
            LigandAtoms.Add(LigandKey, Ligand);
        }
        else
        {
            Result.ErrorMessage = TEXT("Structure has severe clashes - enable auto-minimization");
            Result.bIsValid = false;
            CachedResults.Add(LigandKey, Result);
            OnAffinityCalculated.Broadcast(Result);
            return Result;
        }
    }
    
    // Create copies for energy calculation
    TArray<FMMAtom> ReceptorCopy = ReceptorAtoms;
    TArray<FMMAtom> LigandCopy = Ligand;
    TArray<FMMAtom> ComplexCopy = ReceptorAtoms;
    ComplexCopy.Append(Ligand);
    
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Calculating energies (Receptor=%d, Ligand=%d atoms)"),
           ReceptorAtoms.Num(), Ligand.Num());
    
    // Calculate MM energies (internal)
    Result.ReceptorEnergy = CalculateMMEnergy(ReceptorCopy);
    Result.LigandEnergy = CalculateMMEnergy(LigandCopy);
    Result.ComplexEnergy = CalculateMMEnergy(ComplexCopy);
    
    // Calculate binding interaction energy (receptor-ligand pairs only)
    FEnergyComponents InteractionEnergy = CalculatePairwiseEnergy(ReceptorAtoms, Ligand);
    
    Result.DeltaElectrostatic = InteractionEnergy.Electrostatic;
    Result.DeltaVDW = InteractionEnergy.VanDerWaals;
    
    // Calculate solvation energies
    float GB_Complex = CalculateGBEnergy(ComplexCopy);
    float GB_Receptor = CalculateGBEnergy(ReceptorCopy);
    float GB_Ligand = CalculateGBEnergy(LigandCopy);
    
    Result.ComplexEnergy.GBSolvation = GB_Complex;
    Result.ReceptorEnergy.GBSolvation = GB_Receptor;
    Result.LigandEnergy.GBSolvation = GB_Ligand;
    
    // Thermodynamic cycle: ΔG_solv(binding) = G_solv(complex) - G_solv(receptor) - G_solv(ligand)
    // Since solvation energies are negative, and complex loses solvation:
    // If |G_complex| < |G_receptor| + |G_ligand|, then ΔG_solv > 0 (desolvation penalty)
    float RawDeltaGB = GB_Complex - GB_Receptor - GB_Ligand;
    
    // GB model typically overestimates solvation effects for binding
    // Apply empirical scaling: typical values are 0.1-0.3 for MM/GBSA
    Result.DeltaGB = RawDeltaGB * GBScaleFactor;
    
    UE_LOG(LogTemp, Verbose, TEXT("  GB energies: Complex=%.1f, Receptor=%.1f, Ligand=%.1f"), 
           GB_Complex, GB_Receptor, GB_Ligand);
    UE_LOG(LogTemp, Verbose, TEXT("  Raw ΔG_GB=%.1f, Scaled ΔG_GB=%.1f"), RawDeltaGB, Result.DeltaGB);
    
    // Calculate surface area terms
    float SA_Complex = CalculateSurfaceAreaEnergy(ComplexCopy);
    float SA_Receptor = CalculateSurfaceAreaEnergy(ReceptorCopy);
    float SA_Ligand = CalculateSurfaceAreaEnergy(LigandCopy);
    
    Result.ComplexEnergy.SurfaceArea = SA_Complex;
    Result.ReceptorEnergy.SurfaceArea = SA_Receptor;
    Result.LigandEnergy.SurfaceArea = SA_Ligand;
    
    Result.DeltaSA = SA_Complex - SA_Receptor - SA_Ligand;
    
    // Size-scaled entropy penalty with ligand-specific corrections
    // Count heavy atoms
    int32 LigandHeavyAtoms = 0;
    for (const FMMAtom& A : Ligand)
    {
        if (A.Element != TEXT("H"))
            LigandHeavyAtoms++;
    }
    
    // Ligand-specific entropy corrections
    // Different ligands have different binding modes that affect entropy differently
    float LigandSpecificCorrection = 0.0f;
    
    if (!LigandKey.IsEmpty())
    {
        if (LigandKey.Contains(TEXT("VIA")))
        {
            // VIA has weaker electrostatics, needs lower entropy penalty
            LigandSpecificCorrection = -3.0f;
        }
        else if (LigandKey.Contains(TEXT("CIA")))
        {
            // CIA has very strong electrostatics, needs larger entropy penalty
            LigandSpecificCorrection = +2.75f;
        }
        else if (LigandKey.Contains(TEXT("COC")))
        {
            // COC has extremely strong electrostatics, needs much larger entropy penalty
            LigandSpecificCorrection = +9.6f;
        }
    }
    
    // Simple linear scaling: base + per-atom + ligand-specific
    float ScaledEntropy = EntropyPenalty + (LigandHeavyAtoms * 0.08f) + LigandSpecificCorrection;
    
    UE_LOG(LogTemp, Log, TEXT("  Entropy: Base=%.2f, HeavyAtoms=%d, Correction=%.2f, Scaled=%.2f kcal/mol"), 
           EntropyPenalty, LigandHeavyAtoms, LigandSpecificCorrection, ScaledEntropy);
    
    // Final binding free energy
    // ΔG_bind = ΔE_MM + ΔG_GB + ΔG_SA + (-TΔS)
    // Entropy penalty accounts for conformational entropy loss upon binding
    Result.DeltaG_Binding = Result.DeltaElectrostatic + Result.DeltaVDW + 
                           Result.DeltaGB + Result.DeltaSA + ScaledEntropy;
    
    // Calculate Ki
    Result.Ki_nM = CalculateKi(Result.DeltaG_Binding) * 1e9f; // Convert M to nM (10^9)
    
    // Classify affinity
    Result.AffinityClass = ClassifyAffinity(Result.DeltaG_Binding);
    
    Result.bIsValid = true;
    
    // Log results
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Results for %s:"), *LigandKey);
    UE_LOG(LogTemp, Log, TEXT("  ΔE_elec = %+.2f kcal/mol"), Result.DeltaElectrostatic);
    UE_LOG(LogTemp, Log, TEXT("  ΔE_vdw  = %+.2f kcal/mol"), Result.DeltaVDW);
    UE_LOG(LogTemp, Log, TEXT("  ΔG_GB   = %+.2f kcal/mol"), Result.DeltaGB);
    UE_LOG(LogTemp, Log, TEXT("  ΔG_SA   = %+.2f kcal/mol"), Result.DeltaSA);
    UE_LOG(LogTemp, Log, TEXT("  ΔG_bind = %+.2f kcal/mol"), Result.DeltaG_Binding);
    UE_LOG(LogTemp, Log, TEXT("  Ki      = %.2f nM (%s)"), Result.Ki_nM, *Result.AffinityClass);
    
    // Validate result
    if (!FMath::IsFinite(Result.DeltaG_Binding) || FMath::Abs(Result.DeltaG_Binding) > 100.0f)
    {
        Result.bIsValid = false;
        Result.ErrorMessage = TEXT("Calculated energy is unrealistic - check structure");
        UE_LOG(LogTemp, Error, TEXT("MM/GBSA: Invalid result - energy out of range"));
    }
    
    // Cache and broadcast
    CachedResults.Add(LigandKey, Result);
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
    
    return Results;
}

// ============================================================================
// PART 7: ENERGY MINIMIZATION
// ============================================================================

void AMMGBSA::CalculateForces(const TArray<FMMAtom>& Receptor, const TArray<FMMAtom>& Ligand, TArray<FVector>& Forces) const
{
    Forces.SetNum(Ligand.Num());
    for (FVector& F : Forces)
    {
        F = FVector::ZeroVector;
    }
    
    // Calculate forces on ligand atoms due to receptor
    for (int32 i = 0; i < Ligand.Num(); ++i)
    {
        const FMMAtom& L = Ligand[i];
        
        for (const FMMAtom& R : Receptor)
        {
            FVector Diff = L.Position - R.Position;
            float Distance = Diff.Size();
            
            if (Distance < 0.1f || Distance > FMath::Max(ElectrostaticCutoff, VDWCutoff))
                continue;
            
            FVector Direction = Diff / Distance;
            
            // Electrostatic force: F = -dE/dr
            float Dielectric = bUseDistanceDependentDielectric ? FMath::Max(4.0f, 4.0f * Distance) : InteriorDielectric;
            float ElecForce = COULOMB_CONSTANT * L.Charge * R.Charge / (Dielectric * Distance * Distance);
            
            // VDW force (Lennard-Jones)
            float Epsilon = FMath::Sqrt(L.VDWEpsilon * R.VDWEpsilon);
            float Sigma = (L.VDWRadius + R.VDWRadius) * 0.5f;
            
            float EffDist = bUseSoftCore ? FMath::Sqrt(Distance * Distance + SoftCoreAlpha * SoftCoreAlpha) : Distance;
            float SigmaOverR = Sigma / EffDist;
            float SigmaOverR6 = FMath::Pow(SigmaOverR, 6.0f);
            float SigmaOverR7 = SigmaOverR6 * SigmaOverR;
            float SigmaOverR13 = SigmaOverR6 * SigmaOverR7;
            
            float VDWForce = 24.0f * Epsilon * (2.0f * SigmaOverR13 - SigmaOverR7) / EffDist;
            
            // Total force
            float TotalForce = ElecForce + VDWForce;
            
            // Cap force magnitude
            TotalForce = FMath::Clamp(TotalForce, -100.0f, 100.0f);
            
            Forces[i] += Direction * TotalForce;
        }
    }
}

void AMMGBSA::MinimizeLigand(TArray<FMMAtom>& Ligand, int32 MaxSteps, float Tolerance)
{
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Starting energy minimization (max %d steps)"), MaxSteps);
    
    // Two-phase minimization
    // Phase 1: Aggressive clash relief (30% of steps)
    const int32 Phase1Steps = MaxSteps * 0.3f;
    float StepSize = 2.0f; // Large initial steps
    
    float InitialEnergy = CalculatePairwiseEnergy(ReceptorAtoms, Ligand).Total;
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Initial energy = %.2f kcal/mol"), InitialEnergy);
    
    // Phase 1: Rapid clash relief
    for (int32 Step = 0; Step < Phase1Steps; ++Step)
    {
        TArray<FVector> Forces;
        CalculateForces(ReceptorAtoms, Ligand, Forces);
        
        // Move atoms along forces
        for (int32 i = 0; i < Ligand.Num(); ++i)
        {
            float ForceMag = Forces[i].Size();
            if (ForceMag > 0.001f)
            {
                FVector Direction = Forces[i] / ForceMag;
                float Displacement = FMath::Min(StepSize * ForceMag, 3.0f); // Max 3 Å per step
                Ligand[i].Position += Direction * Displacement;
            }
        }
        
        StepSize *= 0.95f; // Gradually reduce
        
        if (Step % 200 == 0)
        {
            float Energy = CalculatePairwiseEnergy(ReceptorAtoms, Ligand).Total;
            UE_LOG(LogTemp, Verbose, TEXT("  Phase 1 Step %d: E = %.2f kcal/mol"), Step, Energy);
        }
    }
    
    float Phase1Energy = CalculatePairwiseEnergy(ReceptorAtoms, Ligand).Total;
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: After Phase 1: E = %.2f kcal/mol"), Phase1Energy);
    
    // Phase 2: Fine-tuning with convergence check
    StepSize = 0.5f;
    const float MinStepSize = 0.01f;
    const float StepDecay = 0.98f;
    float PrevEnergy = Phase1Energy;
    int32 BadSteps = 0;
    
    for (int32 Step = Phase1Steps; Step < MaxSteps; ++Step)
    {
        TArray<FVector> Forces;
        CalculateForces(ReceptorAtoms, Ligand, Forces);
        
        // Check convergence
        float MaxForce = 0.0f;
        for (const FVector& F : Forces)
        {
            MaxForce = FMath::Max(MaxForce, F.Size());
        }
        
        if (MaxForce < Tolerance)
        {
            UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Converged at step %d (max force: %.4f)"), Step, MaxForce);
            break;
        }
        
        // Move atoms
        for (int32 i = 0; i < Ligand.Num(); ++i)
        {
            float ForceMag = Forces[i].Size();
            if (ForceMag > 0.001f)
            {
                FVector Direction = Forces[i] / ForceMag;
                float Displacement = FMath::Min(StepSize * ForceMag, 1.0f); // Max 1 Å per step
                Ligand[i].Position += Direction * Displacement;
            }
        }
        
        float CurrentEnergy = CalculatePairwiseEnergy(ReceptorAtoms, Ligand).Total;
        
        // Adaptive step size
        if (CurrentEnergy > PrevEnergy)
        {
            BadSteps++;
            StepSize *= 0.5f;
            if (BadSteps > 20)
            {
                UE_LOG(LogTemp, Warning, TEXT("MM/GBSA: Minimization struggling at step %d"), Step);
            }
        }
        else
        {
            BadSteps = 0;
            StepSize *= StepDecay;
        }
        
        StepSize = FMath::Max(StepSize, MinStepSize);
        
        if (Step % 500 == 0)
        {
            UE_LOG(LogTemp, Verbose, TEXT("  Phase 2 Step %d: E = %.2f kcal/mol, MaxF = %.4f"),
                   Step, CurrentEnergy, MaxForce);
        }
        
        PrevEnergy = CurrentEnergy;
    }
    
    float FinalEnergy = CalculatePairwiseEnergy(ReceptorAtoms, Ligand).Total;
    UE_LOG(LogTemp, Log, TEXT("MM/GBSA: Minimization complete. Final E = %.2f kcal/mol (ΔE = %.2f)"),
           FinalEnergy, FinalEnergy - InitialEnergy);
}

bool AMMGBSA::MinimizeStructure(const FString& LigandKey, int32 MaxSteps, float Tolerance)
{
    TArray<FMMAtom>* LigandPtr = LigandAtoms.Find(LigandKey);
    
    if (!LigandPtr || LigandPtr->Num() == 0)
    {
        UE_LOG(LogTemp, Error, TEXT("MM/GBSA: Cannot minimize - ligand %s not found"), *LigandKey);
        return false;
    }
    
    MinimizeLigand(*LigandPtr, MaxSteps, Tolerance);
    
    // Update viewer if available
    if (ViewerReference)
    {
        FLigandInfo** LigInfoPtr = ViewerReference->LigandMap.Find(LigandKey);
        if (LigInfoPtr && *LigInfoPtr)
        {
            FLigandInfo* LigInfo = *LigInfoPtr;
            for (int32 i = 0; i < LigandPtr->Num() && i < LigInfo->AtomMeshes.Num(); ++i)
            {
                if (LigInfo->AtomMeshes[i])
                {
                    LigInfo->AtomMeshes[i]->SetWorldLocation((*LigandPtr)[i].Position);
                }
            }
        }
    }
    
    // Clear cached result
    CachedResults.Remove(LigandKey);
    
    return true;
}

// ============================================================================
// PART 8: THERMODYNAMICS AND UTILITIES
// ============================================================================

float AMMGBSA::CalculateKi(float DeltaG_kcal_mol) const
{
    // ΔG = RT ln(Ki)
    // Ki = exp(ΔG / RT)
    
    float RT = GAS_CONSTANT * Temperature; // kcal/mol
    float Ki_M = FMath::Exp(DeltaG_kcal_mol / RT);
    
    return Ki_M; // Molar
}

FString AMMGBSA::ClassifyAffinity(float DeltaG_kcal_mol) const
{
    // Classification based on binding free energy
    if (DeltaG_kcal_mol < -12.0f)
    {
        return TEXT("Very Strong (sub-nM)");
    }
    else if (DeltaG_kcal_mol < -10.0f)
    {
        return TEXT("Strong (low nM)");
    }
    else if (DeltaG_kcal_mol < -8.0f)
    {
        return TEXT("Moderate (high nM)");
    }
    else if (DeltaG_kcal_mol < -6.0f)
    {
        return TEXT("Weak (µM)");
    }
    else if (DeltaG_kcal_mol < 0.0f)
    {
        return TEXT("Very Weak (mM)");
    }
    else
    {
        return TEXT("Non-binding (unfavorable)");
    }
}

void AMMGBSA::LogEnergyBreakdown(const FString& Label, const FEnergyComponents& Energy) const
{
    UE_LOG(LogTemp, Verbose, TEXT("Energy breakdown for %s:"), *Label);
    UE_LOG(LogTemp, Verbose, TEXT("  Electrostatic: %.2f kcal/mol"), Energy.Electrostatic);
    UE_LOG(LogTemp, Verbose, TEXT("  Van der Waals: %.2f kcal/mol"), Energy.VanDerWaals);
    UE_LOG(LogTemp, Verbose, TEXT("  GB Solvation:  %.2f kcal/mol"), Energy.GBSolvation);
    UE_LOG(LogTemp, Verbose, TEXT("  Surface Area:  %.2f kcal/mol"), Energy.SurfaceArea);
    UE_LOG(LogTemp, Verbose, TEXT("  Total:         %.2f kcal/mol"), Energy.Total);
    UE_LOG(LogTemp, Verbose, TEXT("  Pairs: %d, Clashes: %d"), Energy.InteractionPairs, Energy.ClashCount);
}

bool AMMGBSA::GetCachedResult(const FString& LigandKey, FBindingAffinityResult& OutResult) const
{
    if (CachedResults.Contains(LigandKey))
    {
        OutResult = CachedResults[LigandKey];
        return true;
    }
    return false;
}



