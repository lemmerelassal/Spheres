// MMGBSA.cpp - MM/GBSA Implementation
#include "MMGBSA.h"
#include "PDBViewer.h"
#include "Components/StaticMeshComponent.h"
#include "Kismet/GameplayStatics.h"

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
    
    // Load receptor atoms (from ResidueMap - these are ATOM entries)
    for (const auto& Pair : ViewerReference->ResidueMap)
    {
        const FResidueInfo* ResInfo = Pair.Value;
        if (!ResInfo || !ResInfo->bIsVisible) continue;
        
        // Get atom positions from bond meshes (we only draw bonds for protein)
        // We need to infer atom positions - use a different approach
        // Actually, we should extend PDBViewer to expose actual atom data
        // For now, we'll estimate from the stored parsed data
    }
    
    // Better approach: Re-parse to get atom data
    // We'll need to access the CurrentPDBContent from viewer
    // For this implementation, we'll work with what we have
    
    // Load ligand atoms (from LigandMap - these are HETATM entries or SDF)
    for (const auto& Pair : ViewerReference->LigandMap)
    {
        const FString& Key = Pair.Key;
        FLigandInfo* LigInfo = Pair.Value;
        if (!LigInfo) continue;
        
        TArray<FMMAtom> Atoms;
        
        for (UStaticMeshComponent* Mesh : LigInfo->AtomMeshes)
        {
            if (!Mesh) continue;
            
            FMMAtom Atom;
            Atom.Position = Mesh->GetComponentLocation();
            Atom.Element = TEXT("C"); // Default - would need to read from viewer
            Atom.bIsReceptor = false;
            
            const FElementData* Data = ElementDatabase.Find(Atom.Element);
            if (Data)
            {
                Atom.Radius = Data->VDWRadius * 100.0f; // Convert to UE units
                Atom.Charge = Data->TypicalCharge;
            }
            
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
    
    // Check cache
    if (CachedResults.Contains(LigandKey))
    {
        return CachedResults[LigandKey];
    }
    
    // Find ligand atoms
    TArray<FMMAtom>* LigandPtr = LigandAtoms.Find(LigandKey);
    if (!LigandPtr || LigandPtr->Num() == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("MMGBSA: Ligand %s not found"), *LigandKey);
        return Result;
    }
    
    TArray<FMMAtom>& Ligand = *LigandPtr;
    Result.LigandName = LigandKey;
    
    // Create complex (receptor + ligand)
    TArray<FMMAtom> Complex = ReceptorAtoms;
    Complex.Append(Ligand);
    
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: Calculating for %s (%d atoms)"), 
        *LigandKey, Ligand.Num());
    
    // Calculate MM energies
    Result.EMM_Complex = CalculateMMEnergy(Complex);
    Result.EMM_Receptor = CalculateMMEnergy(ReceptorAtoms);
    Result.EMM_Ligand = CalculateMMEnergy(Ligand);
    Result.DeltaEMM = Result.EMM_Complex - Result.EMM_Receptor - Result.EMM_Ligand;
    
    // Calculate GB solvation energies
    TArray<FMMAtom> ComplexCopy = Complex;
    TArray<FMMAtom> ReceptorCopy = ReceptorAtoms;
    TArray<FMMAtom> LigandCopy = Ligand;
    
    Result.GB_Complex = CalculateGBEnergy(ComplexCopy);
    Result.GB_Receptor = CalculateGBEnergy(ReceptorCopy);
    Result.GB_Ligand = CalculateGBEnergy(LigandCopy);
    Result.DeltaGGB = Result.GB_Complex - Result.GB_Receptor - Result.GB_Ligand;
    
    // Calculate surface area terms
    Result.SA_Complex = CalculateSurfaceAreaEnergy(ComplexCopy);
    Result.SA_Receptor = CalculateSurfaceAreaEnergy(ReceptorCopy);
    Result.SA_Ligand = CalculateSurfaceAreaEnergy(LigandCopy);
    Result.DeltaGSA = Result.SA_Complex - Result.SA_Receptor - Result.SA_Ligand;
    
    // Total binding free energy
    Result.DeltaG_Binding = Result.DeltaEMM + Result.DeltaGGB + Result.DeltaGSA;
    
    // Convert to Ki
    Result.Ki_uM = DeltaGToKi(Result.DeltaG_Binding);
    Result.AffinityClass = ClassifyAffinity(Result.DeltaG_Binding);
    
    UE_LOG(LogTemp, Log, TEXT("MMGBSA: %s | ΔG = %.2f kcal/mol | Ki = %.2f µM | %s"), 
        *LigandKey, Result.DeltaG_Binding, Result.Ki_uM, *Result.AffinityClass);
    
    // Cache result
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

float AMMGBSA::CalculateElectrostaticEnergy(const TArray<FMMAtom>& Atoms) const
{
    // Coulomb's law: E = (1/4πε₀) * Σᵢⱼ (qᵢqⱼ / rᵢⱼ)
    // In vacuum: k = 332.0636 kcal·Å/(mol·e²)
    const float CoulombK = 332.0636f / InteriorDielectric;
    float Energy = 0.0f;
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        for (int32 j = i + 1; j < Atoms.Num(); ++j)
        {
            float Dist = FVector::Dist(Atoms[i].Position, Atoms[j].Position) / 100.0f; // Convert to Å
            if (Dist < 0.5f) continue; // Avoid singularity
            
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
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        for (int32 j = i + 1; j < Atoms.Num(); ++j)
        {
            float Dist = FVector::Dist(Atoms[i].Position, Atoms[j].Position) / 100.0f; // Convert to Å
            if (Dist < 0.5f) continue;
            
            // Lorentz-Berthelot combining rules
            float Epsilon = FMath::Sqrt(
                GetVDWEpsilon(Atoms[i].Element) * 
                GetVDWEpsilon(Atoms[j].Element)
            );
            
            float Sigma = (GetVDWRadius(Atoms[i].Element) + 
                          GetVDWRadius(Atoms[j].Element)) * 0.5f;
            
            float SigmaOverR = Sigma / Dist;
            float SigmaOverR6 = FMath::Pow(SigmaOverR, 6.0f);
            float SigmaOverR12 = SigmaOverR6 * SigmaOverR6;
            
            float E = 4.0f * Epsilon * (SigmaOverR12 - SigmaOverR6);
            Energy += E;
        }
    }
    
    return Energy;
}

// Generalized Born implicit solvation
float AMMGBSA::CalculateGBEnergy(TArray<FMMAtom>& Atoms)
{
    // First calculate Born radii
    CalculateBornRadii(Atoms);
    
    // GB energy: ΔG_GB = -½(1/εᵢₙ - 1/εₒᵤₜ) Σᵢⱼ qᵢqⱼ f_GB
    // f_GB = 1/√(rᵢⱼ² + αᵢαⱼexp(-rᵢⱼ²/4αᵢαⱼ))
    
    const float Factor = -0.5f * 332.0636f * 
        (1.0f / InteriorDielectric - 1.0f / ExteriorDielectric);
    
    float Energy = 0.0f;
    
    // Self energy term
    for (const FMMAtom& Atom : Atoms)
    {
        if (Atom.GBRadius > 0.0f)
        {
            Energy += Factor * Atom.Charge * Atom.Charge / Atom.GBRadius;
        }
    }
    
    // Pairwise terms
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        for (int32 j = i + 1; j < Atoms.Num(); ++j)
        {
            float Dist = FVector::Dist(Atoms[i].Position, Atoms[j].Position) / 100.0f; // Å
            float Ai = Atoms[i].GBRadius;
            float Aj = Atoms[j].GBRadius;
            
            if (Ai <= 0.0f || Aj <= 0.0f) continue;
            
            float AiAj = Ai * Aj;
            float Exp = FMath::Exp(-Dist * Dist / (4.0f * AiAj));
            float fGB = 1.0f / FMath::Sqrt(Dist * Dist + AiAj * Exp);
            
            Energy += 2.0f * Factor * Atoms[i].Charge * Atoms[j].Charge * fGB;
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

float AMMGBSA::CalculateSurfaceAreaEnergy(TArray<FMMAtom>& Atoms)
{
    // ΔG_SA = γ * SASA
    CalculateSASA(Atoms);
    
    float TotalSASA = 0.0f;
    for (const FMMAtom& Atom : Atoms)
    {
        TotalSASA += Atom.SASA;
    }
    
    // Convert surface tension from kcal/(mol·Å²) and SASA from Ų to get energy
    return SurfaceTension * TotalSASA;
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
    
    for (FMMAtom& Atom : Atoms)
    {
        Atom.Charge = EstimateChargeFromElement(Atom.Element);
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

float AMMGBSA::DeltaGToKi(float DeltaG_kcal_mol) const
{
    // ΔG = RT ln(Ki)
    // Ki = exp(ΔG / RT)
    // R = 0.001987 kcal/(mol·K)
    
    const float R = 0.001987f;
    float RT = R * Temperature;
    
    // Ki in M
    float Ki_M = FMath::Exp(DeltaG_kcal_mol / RT);
    
    // Convert to µM
    return Ki_M * 1.0e6f;
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