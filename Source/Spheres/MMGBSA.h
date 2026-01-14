// MMGBSA.h - MM/GBSA Free Energy Calculation for Binding Affinity
// Complete rewrite with improved physics and numerical stability
#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "MMGBSA.generated.h"

class APDBViewer;

// Forward declarations
struct FMMAtom;

// Energy calculation result with detailed diagnostics
USTRUCT(BlueprintType)
struct FEnergyComponents
{
    GENERATED_BODY()

    UPROPERTY(BlueprintReadOnly)
    float Electrostatic = 0.0f;
    
    UPROPERTY(BlueprintReadOnly)
    float VanDerWaals = 0.0f;
    
    UPROPERTY(BlueprintReadOnly)
    float GBSolvation = 0.0f;
    
    UPROPERTY(BlueprintReadOnly)
    float SurfaceArea = 0.0f;
    
    UPROPERTY(BlueprintReadOnly)
    float Total = 0.0f;
    
    // Diagnostic info
    UPROPERTY(BlueprintReadOnly)
    int32 InteractionPairs = 0;
    
    UPROPERTY(BlueprintReadOnly)
    int32 ClashCount = 0;
    
    UPROPERTY(BlueprintReadOnly)
    bool bIsValid = false;
};

// Detailed overlap information for diagnostics
USTRUCT(BlueprintType)
struct FAtomicClash
{
    GENERATED_BODY()

    UPROPERTY(BlueprintReadOnly)
    FString Atom1Label;

    UPROPERTY(BlueprintReadOnly)
    FString Atom2Label;

    UPROPERTY(BlueprintReadOnly)
    FVector Position1;

    UPROPERTY(BlueprintReadOnly)
    FVector Position2;

    UPROPERTY(BlueprintReadOnly)
    float Distance;

    UPROPERTY(BlueprintReadOnly)
    float MinDistance;

    UPROPERTY(BlueprintReadOnly)
    float Severity; // 0-1, where 1 is catastrophic

    FAtomicClash()
        : Position1(FVector::ZeroVector)
        , Position2(FVector::ZeroVector)
        , Distance(0.0f)
        , MinDistance(0.0f)
        , Severity(0.0f)
    {}
};

// Structure quality assessment
USTRUCT(BlueprintType)
struct FStructureQuality
{
    GENERATED_BODY()

    UPROPERTY(BlueprintReadOnly)
    bool bIsAcceptable = false;

    UPROPERTY(BlueprintReadOnly)
    int32 SevereClashes = 0;

    UPROPERTY(BlueprintReadOnly)
    int32 ModerateClashes = 0;

    UPROPERTY(BlueprintReadOnly)
    float WorstClashSeverity = 0.0f;

    UPROPERTY(BlueprintReadOnly)
    FString RecommendedAction;

    UPROPERTY(BlueprintReadOnly)
    TArray<FAtomicClash> TopClashes;
};

// Final binding affinity result
USTRUCT(BlueprintType)
struct FBindingAffinityResult
{
    GENERATED_BODY()
    
    UPROPERTY(BlueprintReadOnly)
    FString LigandName;
    
    UPROPERTY(BlueprintReadOnly)
    FString LigandKey;
    
    // Energy breakdown
    UPROPERTY(BlueprintReadOnly)
    FEnergyComponents ComplexEnergy;
    
    UPROPERTY(BlueprintReadOnly)
    FEnergyComponents ReceptorEnergy;
    
    UPROPERTY(BlueprintReadOnly)
    FEnergyComponents LigandEnergy;
    
    // Delta values (binding contributions)
    UPROPERTY(BlueprintReadOnly)
    float DeltaElectrostatic = 0.0f; // ΔE_elec
    
    UPROPERTY(BlueprintReadOnly)
    float DeltaVDW = 0.0f; // ΔE_vdw
    
    UPROPERTY(BlueprintReadOnly)
    float DeltaGB = 0.0f; // ΔG_solv
    
    UPROPERTY(BlueprintReadOnly)
    float DeltaSA = 0.0f; // ΔG_SA
    
    // Final binding free energy and affinity
    UPROPERTY(BlueprintReadOnly)
    float DeltaG_Binding = 0.0f; // kcal/mol
    
    UPROPERTY(BlueprintReadOnly)
    float Ki_nM = 0.0f; // Dissociation constant in nanomolar
    
    UPROPERTY(BlueprintReadOnly)
    FString AffinityClass; // "Very Strong", "Strong", "Moderate", "Weak"
    
    // Structure quality
    UPROPERTY(BlueprintReadOnly)
    FStructureQuality QualityAssessment;
    
    UPROPERTY(BlueprintReadOnly)
    bool bIsValid = false;
    
    UPROPERTY(BlueprintReadOnly)
    FString ErrorMessage;
    
    // Detailed diagnostics (optional, for debugging)
    UPROPERTY(BlueprintReadOnly)
    TArray<FString> CalculationLog;
};

// Internal atom representation
USTRUCT()
struct FMMAtom
{
    GENERATED_BODY()
    
    FVector Position;
    FString Element;
    FString Label; // e.g., "ARG123:NH1" for diagnostics
    
    float Charge;
    float VDWRadius;
    float VDWEpsilon;
    float BornRadius; // For GB calculation
    float SASA; // Solvent accessible surface area
    
    bool bIsReceptor; // True if part of receptor, false if ligand
    int32 SourceIndex; // Index in original structure
    
    FMMAtom()
        : Position(FVector::ZeroVector)
        , Element(TEXT("C"))
        , Label(TEXT(""))
        , Charge(0.0f)
        , VDWRadius(1.7f)
        , VDWEpsilon(0.086f)
        , BornRadius(0.0f)
        , SASA(0.0f)
        , bIsReceptor(true)
        , SourceIndex(-1)
    {}
};

// Force field parameters for each element
USTRUCT()
struct FForceFieldParams
{
    GENERATED_BODY()
    
    float VDWRadius; // Ångströms
    float VDWEpsilon; // kcal/mol
    float BaseCharge; // Typical partial charge
    float BornRadiusScale; // Scaling factor for GB
    float Electronegativity; // For charge refinement
    
    FForceFieldParams()
        : VDWRadius(1.7f)
        , VDWEpsilon(0.086f)
        , BaseCharge(0.0f)
        , BornRadiusScale(1.0f)
        , Electronegativity(2.5f)
    {}
};

DECLARE_DYNAMIC_MULTICAST_DELEGATE_OneParam(FOnBindingAffinityCalculated, const FBindingAffinityResult&, Result);

UCLASS()
class SPHERES_API AMMGBSA : public AActor
{
    GENERATED_BODY()

public:
    AMMGBSA();
    
    // === Public API ===
    
    // Calculate binding affinity for a specific ligand
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    FBindingAffinityResult CalculateBindingAffinity(const FString& LigandKey);
    
    // Calculate for all visible ligands
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    TArray<FBindingAffinityResult> CalculateAllBindingAffinities();
    
    // Initialize from PDB viewer
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    void InitializeFromViewer(APDBViewer* Viewer);
    
    // Check structure quality without full calculation
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    FStructureQuality AssessStructureQuality(const FString& LigandKey);
    
    // Minimize structure to relieve clashes
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    bool MinimizeStructure(const FString& LigandKey, int32 MaxSteps = 1000, float Tolerance = 0.1f);
    
    // Clear cached results
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    void ClearCache();
    
    // Get cached result for a ligand (for backwards compatibility)
    UFUNCTION(BlueprintPure, Category = "MM/GBSA")
    bool GetCachedResult(const FString& LigandKey, FBindingAffinityResult& OutResult) const;
    
    // Event fired when calculation completes
    UPROPERTY(BlueprintAssignable, Category = "MM/GBSA")
    FOnBindingAffinityCalculated OnAffinityCalculated;

protected:
    virtual void BeginPlay() override;
    
    // === Physics Parameters ===
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Physics")
    float Temperature = 298.15f; // Kelvin
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Physics")
    float InteriorDielectric = 10.0f; // Protein interior (higher = weaker electrostatics)
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Physics")
    float ExteriorDielectric = 78.5f; // Water
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Physics")
    float SaltConcentration = 0.15f; // M (physiological)
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Physics")
    float SurfaceTension = 0.0072f; // kcal/(mol·Ų) for nonpolar solvation
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Physics")
    float SolventProbeRadius = 1.4f; // Ų (water)
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Physics")
    float EntropyPenalty = 5.0f; // Base conformational entropy loss (kcal/mol, +0.08 per heavy atom, +ligand corrections)
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Physics")
    float ChargeScalingFactor = 0.95f; // Scale all charges (0.5-1.0 for unminimized structures)
    
    // === Calculation Settings ===
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Settings")
    float ElectrostaticCutoff = 15.0f; // Ų
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Settings")
    float VDWCutoff = 12.0f; // Ų
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Settings")
    float GBCutoff = 20.0f; // Ų
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Settings")
    float GBScaleFactor = 0.15f; // Empirical scaling for GB term (0.1-0.3 typical)
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Settings")
    bool bUseSoftCore = true; // Soften VDW at short range
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Settings")
    float SoftCoreAlpha = 0.5f; // Ų
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Settings")
    bool bUseDistanceDependentDielectric = true; // More realistic for unminimized structures
    
    // === Quality Control ===
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Quality")
    float SevereClashThreshold = 0.7f; // Fraction of sum of VDW radii
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Quality")
    float ModerateClashThreshold = 0.85f; // Fraction of sum of VDW radii
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Quality")
    int32 MaxAllowedSevereClashes = 2; // Before rejecting structure
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Quality")
    bool bAutoMinimizeOnClash = true;
    
    // === Minimization Settings ===
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Minimization")
    int32 DefaultMinimizationSteps = 5000;
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA|Minimization")
    float MinimizationTolerance = 0.01f; // kcal/(mol·Ų)

private:
    // === Internal State ===
    
    UPROPERTY()
    APDBViewer* ViewerReference;
    
    TArray<FMMAtom> ReceptorAtoms;
    TMap<FString, TArray<FMMAtom>> LigandAtoms;
    TMap<FString, FBindingAffinityResult> CachedResults;
    TMap<FString, FForceFieldParams> ForceField;
    
    // === Core Energy Calculation Functions ===
    
    // Calculate total MM energy for a system
    FEnergyComponents CalculateMMEnergy(const TArray<FMMAtom>& Atoms) const;
    
    // Calculate pairwise interaction energy (receptor-ligand only)
    FEnergyComponents CalculatePairwiseEnergy(const TArray<FMMAtom>& Receptor, const TArray<FMMAtom>& Ligand) const;
    
    // Individual energy terms
    float CalculateElectrostatic(const FMMAtom& A, const FMMAtom& B, float Distance) const;
    float CalculateVDW(const FMMAtom& A, const FMMAtom& B, float Distance) const;
    float CalculateGBPair(const FMMAtom& A, const FMMAtom& B, float Distance) const;
    
    // GB-specific calculations
    void CalculateBornRadii(TArray<FMMAtom>& Atoms) const;
    float CalculateGBEnergy(TArray<FMMAtom>& Atoms) const;
    
    // SASA calculation
    void CalculateSASA(TArray<FMMAtom>& Atoms) const;
    float CalculateSurfaceAreaEnergy(TArray<FMMAtom>& Atoms) const;
    
    // === Structure Quality Assessment ===
    
    FStructureQuality AssessQuality(const TArray<FMMAtom>& Receptor, const TArray<FMMAtom>& Ligand) const;
    TArray<FAtomicClash> DetectClashes(const TArray<FMMAtom>& Receptor, const TArray<FMMAtom>& Ligand) const;
    
    // === Structure Preparation ===
    
    void LoadAtomsFromViewer();
    void AssignForceFieldParameters(TArray<FMMAtom>& Atoms);
    void RefinePartialCharges(TArray<FMMAtom>& Atoms);
    
    // === Energy Minimization ===
    
    void MinimizeLigand(TArray<FMMAtom>& Ligand, int32 MaxSteps, float Tolerance);
    void CalculateForces(const TArray<FMMAtom>& Receptor, const TArray<FMMAtom>& Ligand, TArray<FVector>& Forces) const;
    
    // === Thermodynamics ===
    
    float CalculateKi(float DeltaG_kcal_mol) const;
    FString ClassifyAffinity(float DeltaG_kcal_mol) const;
    
    // === Utilities ===
    
    void InitializeForceField();
    FForceFieldParams GetForceFieldParams(const FString& Element) const;
    void LogEnergyBreakdown(const FString& Label, const FEnergyComponents& Energy) const;
    
    UFUNCTION()
    void OnViewerLigandsLoaded();
};
