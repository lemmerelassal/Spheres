// MMGBSA.h - MM/GBSA Free Energy Calculation for Binding Affinity
#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "MMGBSA.generated.h"

class APDBViewer;
class AMolecularDynamics;

USTRUCT(BlueprintType)
struct FMMOverlapInfo
{
    GENERATED_BODY()

    UPROPERTY(BlueprintReadOnly)
    FString ReceptorSourceKey;

    UPROPERTY(BlueprintReadOnly)
    int32 ReceptorAtomIndex;

    UPROPERTY(BlueprintReadOnly)
    FString ReceptorElement;

    UPROPERTY(BlueprintReadOnly)
    FVector ReceptorPosition;

    UPROPERTY(BlueprintReadOnly)
    int32 LigandAtomIndex;

    UPROPERTY(BlueprintReadOnly)
    FString LigandElement;

    UPROPERTY(BlueprintReadOnly)
    FVector LigandPosition;

    UPROPERTY(BlueprintReadOnly)
    float DistanceA;

    UPROPERTY(BlueprintReadOnly)
    float SumRadiiA;

    UPROPERTY(BlueprintReadOnly)
    float PairVDW;

    FMMOverlapInfo()
        : ReceptorSourceKey(TEXT(""))
        , ReceptorAtomIndex(-1)
        , ReceptorElement(TEXT(""))
        , ReceptorPosition(FVector::ZeroVector)
        , LigandAtomIndex(-1)
        , LigandElement(TEXT(""))
        , LigandPosition(FVector::ZeroVector)
        , DistanceA(0.0f)
        , SumRadiiA(0.0f)
        , PairVDW(0.0f)
    {}
};

// Result of MM/GBSA calculation for a single ligand
USTRUCT(BlueprintType)
struct FBindingAffinityResult
{
    GENERATED_BODY()
    
    UPROPERTY(BlueprintReadOnly)
    FString LigandName;
    
    UPROPERTY(BlueprintReadOnly)
    FString LigandKey;
    
    // Energy components (kcal/mol)
    UPROPERTY(BlueprintReadOnly)
    float EMM_Complex;  // MM energy of complex
    
    UPROPERTY(BlueprintReadOnly)
    float EMM_Receptor; // MM energy of receptor alone
    
    UPROPERTY(BlueprintReadOnly)
    float EMM_Ligand;   // MM energy of ligand alone
    
    UPROPERTY(BlueprintReadOnly)
    float DeltaEMM;     // EMM_Complex - EMM_Receptor - EMM_Ligand
    
    UPROPERTY(BlueprintReadOnly)
    float GB_Complex;   // GB solvation of complex
    
    UPROPERTY(BlueprintReadOnly)
    float GB_Receptor;  // GB solvation of receptor
    
    UPROPERTY(BlueprintReadOnly)
    float GB_Ligand;    // GB solvation of ligand
    
    UPROPERTY(BlueprintReadOnly)
    float DeltaGGB;     // GB_Complex - GB_Receptor - GB_Ligand
    
    UPROPERTY(BlueprintReadOnly)
    float SA_Complex;   // Surface area term for complex
    
    UPROPERTY(BlueprintReadOnly)
    float SA_Receptor;  // Surface area term for receptor
    
    UPROPERTY(BlueprintReadOnly)
    float SA_Ligand;    // Surface area term for ligand
    
    UPROPERTY(BlueprintReadOnly)
    float DeltaGSA;     // SA_Complex - SA_Receptor - SA_Ligand
    
    // Final binding free energy
    UPROPERTY(BlueprintReadOnly)
    float DeltaG_Binding; // DeltaEMM + DeltaGGB + DeltaGSA
    
    // Dissociation constant (Ki) in micromolar
    UPROPERTY(BlueprintReadOnly)
    float Ki_uM;
    
    // Binding affinity category
    UPROPERTY(BlueprintReadOnly)
    FString AffinityClass; // "Very Strong", "Strong", "Moderate", "Weak"

    // Per-atom overlap diagnostics
    UPROPERTY(BlueprintReadOnly)
    int32 OverlapCount;

    UPROPERTY(BlueprintReadOnly)
    bool bHasSevereOverlap;

    // Detailed overlap info
    UPROPERTY(BlueprintReadOnly)
    TArray<struct FMMOverlapInfo> Overlaps;

    // Whether this result is valid (calculated)
    UPROPERTY(BlueprintReadOnly)
    bool bIsValid;
    
    FBindingAffinityResult()
        : EMM_Complex(0.0f)
        , EMM_Receptor(0.0f)
        , EMM_Ligand(0.0f)
        , DeltaEMM(0.0f)
        , GB_Complex(0.0f)
        , GB_Receptor(0.0f)
        , GB_Ligand(0.0f)
        , DeltaGGB(0.0f)
        , SA_Complex(0.0f)
        , SA_Receptor(0.0f)
        , SA_Ligand(0.0f)
        , DeltaGSA(0.0f)
        , DeltaG_Binding(0.0f)
        , Ki_uM(0.0f)
        , AffinityClass(TEXT("Unknown"))
        , OverlapCount(0)
        , bHasSevereOverlap(false)
        , bIsValid(false)
    {}
};

// Atom information for energy calculations
USTRUCT()
struct FMMAtom
{
    GENERATED_BODY()
    
    FVector Position;
    FString Element;
    float Charge;
    float Radius;      // Van der Waals radius
    float GBRadius;    // Born radius for GB
    float SASA;        // Solvent accessible surface area
    bool bIsReceptor;  // True if part of receptor, false if ligand
    
    // Source tracking for diagnostics
    FString SourceKey; // Residue key (for receptor) or Ligand key (for ligand)
    int32 SourceIndex; // Atom index within the source

    FMMAtom()
        : Position(FVector::ZeroVector)
        , Element(TEXT("C"))
        , Charge(0.0f)
        , Radius(170.0f)
        , GBRadius(0.0f)
        , SASA(0.0f)
        , bIsReceptor(true)
        , SourceIndex(-1)
    {}
};

DECLARE_DYNAMIC_MULTICAST_DELEGATE_OneParam(FOnBindingAffinityCalculated, const FBindingAffinityResult&, Result);

UCLASS()
class SPHERES_API AMMGBSA : public AActor
{
    GENERATED_BODY()

public:
    AMMGBSA();
    
    // Calculate binding affinity for a specific ligand
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    FBindingAffinityResult CalculateBindingAffinity(const FString& LigandKey);
    
    // Calculate binding affinities for all visible ligands
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    TArray<FBindingAffinityResult> CalculateAllBindingAffinities();
    
    // Initialize from PDB viewer
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    void InitializeFromViewer(APDBViewer* Viewer);
    
    // Get cached result for a ligand
    UFUNCTION(BlueprintPure, Category = "MM/GBSA")
    bool GetCachedResult(const FString& LigandKey, FBindingAffinityResult& OutResult) const;
    
    // Set temperature for free energy calculation
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    void SetTemperature(float TempKelvin) { Temperature = TempKelvin; }
    
    // Event fired when calculation completes
    UPROPERTY(BlueprintAssignable, Category = "MM/GBSA")
    FOnBindingAffinityCalculated OnAffinityCalculated;

    // Run a parameter sweep over GBScale and ChargeScaling for a set of ligands. Results are logged and written to Saved/MMGBSA_Sweep_<timestamp>.csv
    UFUNCTION(BlueprintCallable, Category = "MM/GBSA")
    void RunParameterSweep(const TArray<float>& GBScales, const TArray<float>& ChargeScales, const TArray<FString>& LigandKeys);


protected:
    virtual void BeginPlay() override;
    
    // MM/GBSA parameters
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float Temperature = 298.15f; // Kelvin (25°C)
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float InteriorDielectric = 2.0f; // Protein interior (increased from 1.0 to more realistic 2.0)
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float ExteriorDielectric = 78.5f; // Water
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float SurfaceTension = 0.0072f; // kcal/(mol·Ų)
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float SolventProbeRadius = 1.4f; // Å (water probe)
    
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float SaltConcentration = 0.15f; // M (150 mM, physiological)

    // Soft-core van der Waals options: when enabled, pairwise VDW uses a softened distance to avoid huge repulsive spikes
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    bool bUseSoftCoreVDW = true;

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float SoftCoreAlpha = 0.5f; // Å - soft-core smoothing parameter (r_sc = sqrt(r^2 + alpha^2))

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float MaxPairVDW = 100.0f; // kcal/mol - cap for pairwise VDW to avoid numeric explosion

    // Receptor-ligand specific cutoff (Å) for pairwise contributions
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float ReceptorLigandCutoffA = 6.0f; // Å - only consider RL pairs within this distance

    // GB-specific parameters to avoid large solvation spikes
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float ReceptorLigandGBCutoffA = 8.0f; // Å - only consider GB pair contributions within this distance

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float MaxGBPerAtom = 200.0f; // kcal/mol - cap for per-atom GB self-term

    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float MaxGBPair = 200.0f; // kcal/mol - cap for GB pair contributions

    // Charge scaling to reduce magnitude of assigned typical charges (0.0-1.0)
    // Lower values (0.3-0.5) are more realistic since typical charges are approximate
    // and real protein partial charges are usually smaller than element-typical values
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float ChargeScaling = 0.4f; // scale factor applied to typical atomic charges (reduced for realistic values)

    // Minimum allowed GB radius (Å) to avoid tiny radii causing huge self-energies
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float MinGBRadius = 0.8f; // Å

    // Clamp for ΔGGB magnitude to protect against numeric explosion (kcal/mol)
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float MaxDeltaGB = 500.0f;

    // GB scaling factor applied to computed ΔGGB (useful as an empirical correction)
    // Lower values (0.02-0.1) are typically needed to get realistic binding energies
    // because pairwise GB only captures receptor-ligand interactions, not the full
    // solvation change (which includes desolvation penalties that offset the interactions)
    UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "MM/GBSA Parameters")
    float GBScale = 0.03f; // multiply ΔGGB by this factor before clamping (empirical scaling for realistic values)

private:
    UPROPERTY()
    APDBViewer* ViewerReference;
    
    // Cached atom data
    TArray<FMMAtom> ReceptorAtoms;
    TMap<FString, TArray<FMMAtom>> LigandAtoms;
    
    // Cached results
    TMap<FString, FBindingAffinityResult> CachedResults;
    
    // Energy calculation methods
    float CalculateMMEnergy(const TArray<FMMAtom>& Atoms) const;
    float CalculateMMEnergyPair(const TArray<FMMAtom>& Atoms1, const TArray<FMMAtom>& Atoms2) const;
    float CalculateElectrostaticEnergy(const TArray<FMMAtom>& Atoms) const;
    float CalculateVanDerWaalsEnergy(const TArray<FMMAtom>& Atoms) const;
    
    // Generalized Born solvation
    float CalculateGBEnergy(const TArray<FMMAtom>& Atoms);
    float CalculateGBEnergyPair(const TArray<FMMAtom>& Atoms1, const TArray<FMMAtom>& Atoms2); // NEW: pairwise GB only
    void CalculateBornRadii(TArray<FMMAtom>& Atoms);
    float GetEffectiveBornRadius(const FMMAtom& Atom, const TArray<FMMAtom>& AllAtoms);
    
    // Surface area term
    float CalculateSurfaceAreaEnergy(const TArray<FMMAtom>& Atoms);
    void CalculateSASA(TArray<FMMAtom>& Atoms);
    
    // Helper functions
    void LoadAtomsFromViewer();
    void AssignPartialCharges(TArray<FMMAtom>& Atoms);

    // Called when viewer loads/refreshes ligands
    UFUNCTION()
    void OnViewerLigandsLoaded();
    float EstimateChargeFromElement(const FString& Element) const;
    float GetVDWRadius(const FString& Element) const;
    float GetVDWEpsilon(const FString& Element) const;
    
    // Convert ΔG to Ki
    float DeltaGToKi(float DeltaG_kcal_mol) const;
    FString ClassifyAffinity(float DeltaG_kcal_mol) const;
    
    // Element database
    struct FElementData
    {
        float VDWRadius;
        float VDWEpsilon;
        float TypicalCharge;
        float GBRadiusScale;
    };
    TMap<FString, FElementData> ElementDatabase;
    void InitializeElementDatabase();
};