// FEPCalculator.cpp - Free Energy Perturbation Implementation
// Compatible with Unreal Engine 5.6

#include "FEPCalculator.h"
#include "FEPControlWidget.h"
#include "PDBViewer.h"
#include "Components/StaticMeshComponent.h"
#include "Materials/MaterialInstanceDynamic.h"
#include "Engine/StaticMesh.h"
#include "UObject/ConstructorHelpers.h"
#include "TimerManager.h"
#include "Math/UnrealMathUtility.h"
#include "Kismet/GameplayStatics.h"
#include "Misc/FileHelper.h"
#include "Misc/Paths.h"
#include "HAL/PlatformFilemanager.h"

// Physical constants
namespace FEPConstants
{
    constexpr float BOLTZMANN = 0.001987204f; // kcal/(mol·K)
    constexpr float COULOMB = 332.0636f; // e²·Å·kcal/mol (for electrostatics)
    constexpr float AVOGADRO = 6.02214076e23f;
    constexpr float GAS_CONSTANT = 1.987204e-3f; // kcal/(mol·K)
}

AFEPCalculator::AFEPCalculator()
{
    PrimaryActorTick.bCanEverTick = true;
    RootComponent = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
    
    bIsCalculating = false;
    CurrentProgress = 0.0f;
    CurrentLambdaIndex = 0;
    CurrentStep = 0;
    CurrentQueueIndex = 0;
    PDBViewer = nullptr;
}

void AFEPCalculator::BeginPlay()
{
    Super::BeginPlay();
    
    // Auto-initialize with PDBViewer
    InitializeWithPDBViewer();
    
    // Auto-spawn UI widget if enabled
    if (bAutoSpawnUI && ControlWidgetClass)
    {
        APlayerController* PC = GetWorld()->GetFirstPlayerController();
        if (PC)
        {
            ControlWidget = CreateWidget<UFEPControlWidget>(PC, ControlWidgetClass);
            if (ControlWidget)
            {
                ControlWidget->SetFEPCalculator(this);
                ControlWidget->SetPDBViewer(PDBViewer);
                ControlWidget->AddToViewport();
                
                UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Auto-spawned UI widget"));
            }
        }
    }
    
    // Auto-calculate if enabled
    if (bAutoCalculateOnBeginPlay)
    {
        // Delay slightly to ensure everything is initialized
        FTimerHandle DelayHandle;
        GetWorld()->GetTimerManager().SetTimer(DelayHandle, [this]()
        {
            if (bAutoCalculateVisibleOnly)
            {
                UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Starting auto-calculation for visible ligands"));
                CalculateVisibleLigands();
            }
            else if (!AutoCalculateLigandKey.IsEmpty())
            {
                UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Starting auto-calculation for %s"), *AutoCalculateLigandKey);
                CalculateBindingFreeEnergy(AutoCalculateLigandKey);
            }
            else
            {
                UE_LOG(LogTemp, Warning, TEXT("FEPCalculator: Auto-calculate enabled but no ligand specified"));
            }
        }, 0.5f, false);
    }
}

void AFEPCalculator::InitializeWithPDBViewer()
{
    // Find PDBViewer in the level
    TArray<AActor*> FoundActors;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), APDBViewer::StaticClass(), FoundActors);
    
    if (FoundActors.Num() > 0)
    {
        PDBViewer = Cast<APDBViewer>(FoundActors[0]);
        if (PDBViewer)
        {
            UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Found and attached to PDBViewer"));
        }
    }
    else
    {
        UE_LOG(LogTemp, Warning, TEXT("FEPCalculator: No PDBViewer found in level. Place a PDBViewer actor first."));
    }
}

void AFEPCalculator::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);
}

void AFEPCalculator::CalculateBindingFreeEnergy(const FString& LigandKey)
{
    if (!PDBViewer)
    {
        UE_LOG(LogTemp, Error, TEXT("FEPCalculator: PDBViewer not initialized. Call InitializeWithPDBViewer() first."));
        return;
    }
    
    if (bIsCalculating)
    {
        UE_LOG(LogTemp, Warning, TEXT("FEPCalculator: Calculation already in progress"));
        return;
    }
    
    UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Starting binding free energy calculation for ligand: %s"), *LigandKey);
    
    // Initialize
    bIsCalculating = true;
    CurrentProgress = 0.0f;
    CurrentLambdaIndex = 0;
    CurrentStep = 0;
    CollectedEnergies.Empty();
    CollecteddHdL.Empty();
    
    // Initialize system
    InitializeSystem(LigandKey);
    
    // Setup lambda windows
    LastResult.LambdaWindows.Empty();
    for (int32 i = 0; i < Parameters.NumLambdaWindows; ++i)
    {
        FFEPLambdaWindow Window;
        Window.Lambda = float(i) / float(Parameters.NumLambdaWindows - 1);
        LastResult.LambdaWindows.Add(Window);
    }
    
    // Start calculation with timer (async-like processing)
    GetWorld()->GetTimerManager().SetTimer(
        CalculationTimerHandle,
        [this, LigandKey]()
        {
            if (CurrentLambdaIndex < Parameters.NumLambdaWindows)
            {
                float Lambda = LastResult.LambdaWindows[CurrentLambdaIndex].Lambda;
                RunLambdaWindow(Lambda);
            }
            else
            {
                // All windows complete - finalize
                GetWorld()->GetTimerManager().ClearTimer(CalculationTimerHandle);
                
                // Calculate free energy
                LastResult.DeltaG = IntegrateFreeEnergy(LastResult.LambdaWindows);
                LastResult.BindingAffinity = CalculateBindingAffinityFromDeltaG(LastResult.DeltaG);
                LastResult.bCalculationSuccessful = true;
                
                bIsCalculating = false;
                CurrentProgress = 1.0f;
                
                // Store result
                AllResults.Add(LigandKey, LastResult);
                
                UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Calculation complete for %s!"), *LigandKey);
                UE_LOG(LogTemp, Log, TEXT("  ΔG_bind = %.2f ± %.2f kcal/mol"), 
                       LastResult.DeltaG, LastResult.StandardError);
                UE_LOG(LogTemp, Log, TEXT("  K_d = %.2f nM"), LastResult.BindingAffinity);
                
                // Auto-export if enabled
                if (bAutoExportResults)
                {
                    ExportResult(LigandKey, LastResult);
                }
                
                OnFEPComplete.Broadcast(LastResult);
                
                // Process next in queue if any
                ProcessNextInQueue();
            }
        },
        0.016f, // ~60 FPS
        true
    );
}

void AFEPCalculator::CalculateAllLigands()
{
    if (!PDBViewer)
    {
        UE_LOG(LogTemp, Error, TEXT("FEPCalculator: PDBViewer not initialized"));
        return;
    }
    
    if (bIsCalculating)
    {
        UE_LOG(LogTemp, Warning, TEXT("FEPCalculator: Calculation already in progress"));
        return;
    }
    
    // Get all ligand keys
    TArray<FString> LigandKeys;
    PDBViewer->LigandMap.GetKeys(LigandKeys);
    
    if (LigandKeys.Num() == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("FEPCalculator: No ligands found in PDBViewer"));
        return;
    }
    
    UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Starting batch calculation for %d ligands"), LigandKeys.Num());
    
    // Setup queue
    CalculationQueue = LigandKeys;
    CurrentQueueIndex = 0;
    
    // Start first calculation
    if (CalculationQueue.Num() > 0)
    {
        CalculateBindingFreeEnergy(CalculationQueue[0]);
    }
}

void AFEPCalculator::CalculateVisibleLigands()
{
    if (!PDBViewer)
    {
        UE_LOG(LogTemp, Error, TEXT("FEPCalculator: PDBViewer not initialized"));
        return;
    }
    
    if (bIsCalculating)
    {
        UE_LOG(LogTemp, Warning, TEXT("FEPCalculator: Calculation already in progress"));
        return;
    }
    
    // Get visible ligands only
    TArray<FString> VisibleLigandKeys;
    for (auto& Pair : PDBViewer->LigandMap)
    {
        if (Pair.Value && Pair.Value->bIsVisible)
        {
            VisibleLigandKeys.Add(Pair.Key);
        }
    }
    
    if (VisibleLigandKeys.Num() == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("FEPCalculator: No visible ligands found"));
        UE_LOG(LogTemp, Log, TEXT("  Make sure to show at least one ligand in PDBViewer before calculating"));
        return;
    }
    
    UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Starting calculation for %d visible ligand(s)"), VisibleLigandKeys.Num());
    
    // Setup queue
    CalculationQueue = VisibleLigandKeys;
    CurrentQueueIndex = 0;
    
    // Start first calculation
    if (CalculationQueue.Num() > 0)
    {
        CalculateBindingFreeEnergy(CalculationQueue[0]);
    }
}

void AFEPCalculator::ProcessNextInQueue()
{
    CurrentQueueIndex++;
    
    if (CurrentQueueIndex < CalculationQueue.Num())
    {
        UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Processing queue item %d/%d"), 
               CurrentQueueIndex + 1, CalculationQueue.Num());
        
        // Small delay between calculations
        FTimerHandle DelayHandle;
        GetWorld()->GetTimerManager().SetTimer(DelayHandle, [this]()
        {
            CalculateBindingFreeEnergy(CalculationQueue[CurrentQueueIndex]);
        }, 0.5f, false);
    }
    else
    {
        // Queue complete
        if (CalculationQueue.Num() > 0)
        {
            UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Batch calculation complete!"));
            UE_LOG(LogTemp, Log, TEXT("  Calculated %d ligands"), AllResults.Num());
            
            // Sort and display results
            TArray<FString> Keys;
            AllResults.GetKeys(Keys);
            
            Keys.Sort([this](const FString& A, const FString& B)
            {
                return AllResults[A].DeltaG < AllResults[B].DeltaG;
            });
            
            UE_LOG(LogTemp, Log, TEXT("  Ranking (best to worst):"));
            int32 Rank = 1;
            for (const FString& Key : Keys)
            {
                const FFEPResult& Result = AllResults[Key];
                UE_LOG(LogTemp, Log, TEXT("    %d. %s: ΔG = %.2f kcal/mol, Kd = %.2f nM"), 
                       Rank++, *Key, Result.DeltaG, Result.BindingAffinity);
            }
        }
        
        CalculationQueue.Empty();
        CurrentQueueIndex = 0;
    }
}

void AFEPCalculator::CalculateRelativeBindingFreeEnergy(const FString& LigandKeyA, 
                                                        const FString& LigandKeyB)
{
    if (!PDBViewer)
    {
        UE_LOG(LogTemp, Error, TEXT("FEPCalculator: PDBViewer not initialized"));
        return;
    }
    
    UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Relative binding free energy calculation"));
    UE_LOG(LogTemp, Log, TEXT("  This will calculate ΔΔG = ΔG_B - ΔG_A"));
    
    // Add both to queue
    CalculationQueue.Empty();
    CalculationQueue.Add(LigandKeyA);
    CalculationQueue.Add(LigandKeyB);
    CurrentQueueIndex = 0;
    
    // Calculate first ligand
    CalculateBindingFreeEnergy(LigandKeyA);
    
    // After both complete, the results will be in AllResults
    // User can access via GetAllResults() or check logs for ranking
}

void AFEPCalculator::StopCalculation()
{
    if (bIsCalculating)
    {
        GetWorld()->GetTimerManager().ClearTimer(CalculationTimerHandle);
        bIsCalculating = false;
        
        LastResult.bCalculationSuccessful = false;
        LastResult.ErrorMessage = TEXT("Calculation stopped by user");
        
        UE_LOG(LogTemp, Warning, TEXT("FEPCalculator: Calculation stopped"));
    }
}

void AFEPCalculator::InitializeSystem(const FString& LigandKey)
{
    if (!PDBViewer)
    {
        UE_LOG(LogTemp, Error, TEXT("FEPCalculator: PDBViewer not set"));
        return;
    }
    
    CurrentState.Atoms.Empty();
    InitialProteinPositions.Empty();
    
    // Get ligand info
    FLigandInfo* Ligand = nullptr;
    if (PDBViewer->LigandMap.Contains(LigandKey))
    {
        Ligand = PDBViewer->LigandMap[LigandKey];
    }
    
    if (!Ligand)
    {
        UE_LOG(LogTemp, Error, TEXT("FEPCalculator: Ligand not found: %s"), *LigandKey);
        return;
    }
    
    UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Initializing system with %d ligand atoms"), 
           Ligand->AtomPositions.Num());
    
    // Add ligand atoms
    for (int32 i = 0; i < Ligand->AtomPositions.Num(); ++i)
    {
        FAtomState Atom;
        Atom.Position = Ligand->AtomPositions[i];
        Atom.Element = Ligand->AtomElements.IsValidIndex(i) ? Ligand->AtomElements[i] : TEXT("C");
        Atom.bIsLigandAtom = true;
        
        SetAtomParameters(Atom, Atom.Element);
        
        // Initialize velocity from Maxwell-Boltzmann distribution
        float StdDev = FMath::Sqrt(FEPConstants::BOLTZMANN * Parameters.Temperature / Atom.Mass);
        Atom.Velocity.X = FMath::FRandRange(-StdDev, StdDev);
        Atom.Velocity.Y = FMath::FRandRange(-StdDev, StdDev);
        Atom.Velocity.Z = FMath::FRandRange(-StdDev, StdDev);
        
        CurrentState.Atoms.Add(Atom);
    }
    
    // Add nearby protein atoms (within cutoff of ligand)
    FVector LigandCenter = FVector::ZeroVector;
    for (const FVector& Pos : Ligand->AtomPositions)
    {
        LigandCenter += Pos;
    }
    LigandCenter /= Ligand->AtomPositions.Num();
    
    int32 ProteinAtomCount = 0;
    for (auto& ResPair : PDBViewer->ResidueMap)
    {
        FResidueInfo* Residue = ResPair.Value;
        if (!Residue) continue;
        
        for (int32 i = 0; i < Residue->AtomPositions.Num(); ++i)
        {
            FVector AtomPos = Residue->AtomPositions[i];
            float Distance = FVector::Dist(AtomPos, LigandCenter);
            
            // Include atoms within cutoff distance
            if (Distance < Parameters.CutoffDistance * 2.0f)
            {
                FAtomState Atom;
                Atom.Position = AtomPos;
                Atom.Element = Residue->AtomElements.IsValidIndex(i) ? 
                               Residue->AtomElements[i] : TEXT("C");
                Atom.bIsLigandAtom = false;
                
                SetAtomParameters(Atom, Atom.Element);
                
                // Protein atoms start at rest or with minimal velocity
                Atom.Velocity = FVector::ZeroVector;
                
                CurrentState.Atoms.Add(Atom);
                InitialProteinPositions.Add(AtomPos);
                ProteinAtomCount++;
            }
        }
    }
    
    UE_LOG(LogTemp, Log, TEXT("FEPCalculator: System initialized with %d total atoms (%d protein)"), 
           CurrentState.Atoms.Num(), ProteinAtomCount);
    
    CurrentState.Temperature = Parameters.Temperature;
    CurrentState.Lambda = 0.0f;
}

void AFEPCalculator::ExportResult(const FString& LigandKey, const FFEPResult& Result)
{
    if (!Result.bCalculationSuccessful)
    {
        return;
    }
    
    // Create output directory
    FString OutputDir = FPaths::ProjectSavedDir() + TEXT("FEP/");
    FString Timestamp = FDateTime::Now().ToString(TEXT("%Y%m%d_%H%M%S"));
    FString Filename = FString::Printf(TEXT("FEP_%s_%s.txt"), *LigandKey, *Timestamp);
    FString FilePath = OutputDir + Filename;
    
    // Create directory if it doesn't exist
    IPlatformFile& PlatformFile = FPlatformFileManager::Get().GetPlatformFile();
    if (!PlatformFile.DirectoryExists(*OutputDir))
    {
        PlatformFile.CreateDirectory(*OutputDir);
    }
    
    // Build output text
    FString Output;
    Output += TEXT("===========================================\n");
    Output += TEXT("Free Energy Perturbation (FEP) Results\n");
    Output += TEXT("===========================================\n\n");
    Output += FString::Printf(TEXT("Ligand: %s\n"), *LigandKey);
    Output += FString::Printf(TEXT("Timestamp: %s\n\n"), *FDateTime::Now().ToString());
    
    Output += TEXT("BINDING FREE ENERGY:\n");
    Output += FString::Printf(TEXT("  ΔG_bind = %.3f ± %.3f kcal/mol\n"), 
                             Result.DeltaG, Result.StandardError);
    Output += FString::Printf(TEXT("  K_d = %.2f nM\n\n"), Result.BindingAffinity);
    
    Output += TEXT("ENERGY COMPONENTS:\n");
    Output += FString::Printf(TEXT("  Electrostatic: %.3f kcal/mol\n"), 
                             Result.ElectrostaticContribution);
    Output += FString::Printf(TEXT("  Van der Waals: %.3f kcal/mol\n"), 
                             Result.VdWContribution);
    Output += FString::Printf(TEXT("  Solvation: %.3f kcal/mol\n\n"), 
                             Result.SolvationContribution);
    
    Output += TEXT("LAMBDA WINDOWS:\n");
    Output += TEXT("  λ        <E>         <dH/dλ>      ±Error     Samples\n");
    Output += TEXT("  -------------------------------------------------------\n");
    
    for (const FFEPLambdaWindow& Window : Result.LambdaWindows)
    {
        Output += FString::Printf(TEXT("  %.3f    %8.2f    %8.2f    %6.2f    %6d\n"),
                                 Window.Lambda,
                                 Window.Energy,
                                 Window.dHdLambda,
                                 Window.StandardError,
                                 Window.SampleCount);
    }
    
    Output += TEXT("\n===========================================\n");
    
    // Write to file
    if (FFileHelper::SaveStringToFile(Output, *FilePath))
    {
        UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Results exported to %s"), *FilePath);
    }
    else
    {
        UE_LOG(LogTemp, Error, TEXT("FEPCalculator: Failed to export results"));
    }
}

void AFEPCalculator::RunLambdaWindow(float Lambda)
{
    if (bVerboseLogging)
    {
        UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Running lambda = %.3f"), Lambda);
    }
    
    CurrentState.Lambda = Lambda;
    CollectedEnergies.Empty();
    CollecteddHdL.Empty();
    
    // Equilibration phase
    for (int32 Step = 0; Step < Parameters.EquilibrationSteps; ++Step)
    {
        PerformMDStep(Lambda);
        
        if (Step % 1000 == 0)
        {
            float Progress = float(CurrentLambdaIndex) / float(Parameters.NumLambdaWindows);
            Progress += (float(Step) / float(Parameters.EquilibrationSteps)) / float(Parameters.NumLambdaWindows);
            CurrentProgress = Progress * 0.5f; // Equilibration is first 50%
            OnFEPProgress.Broadcast(CurrentProgress * 100.0f);
        }
    }
    
    // Production phase - collect data
    for (int32 Step = 0; Step < Parameters.ProductionSteps; ++Step)
    {
        PerformMDStep(Lambda);
        
        if (Step % Parameters.SamplingInterval == 0)
        {
            float Energy = CalculateTotalEnergy(Lambda);
            float dHdL = CalculatedHdLambda(Lambda);
            
            CollectedEnergies.Add(Energy);
            CollecteddHdL.Add(dHdL);
        }
        
        if (Step % 1000 == 0)
        {
            float Progress = float(CurrentLambdaIndex) / float(Parameters.NumLambdaWindows);
            Progress += (float(Step) / float(Parameters.ProductionSteps)) / float(Parameters.NumLambdaWindows);
            CurrentProgress = 0.5f + Progress * 0.5f; // Production is second 50%
            OnFEPProgress.Broadcast(CurrentProgress * 100.0f);
        }
    }
    
    // Store results for this lambda window
    FFEPLambdaWindow& Window = LastResult.LambdaWindows[CurrentLambdaIndex];
    
    // Calculate average energy
    float SumEnergy = 0.0f;
    for (float E : CollectedEnergies)
    {
        SumEnergy += E;
    }
    Window.Energy = SumEnergy / CollectedEnergies.Num();
    
    // Calculate average dH/dλ
    float SumdHdL = 0.0f;
    for (float dH : CollecteddHdL)
    {
        SumdHdL += dH;
    }
    Window.dHdLambda = SumdHdL / CollecteddHdL.Num();
    Window.SampleCount = CollecteddHdL.Num();
    
    // Calculate standard error
    Window.StandardError = CalculateBlockAverageError(CollecteddHdL);
    
    if (bVerboseLogging)
    {
        UE_LOG(LogTemp, Log, TEXT("  <E> = %.2f, <dH/dλ> = %.2f ± %.2f"), 
               Window.Energy, Window.dHdLambda, Window.StandardError);
    }
    
    CurrentLambdaIndex++;
}

void AFEPCalculator::PerformMDStep(float Lambda)
{
    // Velocity Verlet integration
    
    // Calculate forces
    CalculateForces(Lambda);
    
    // Update positions and velocities
    UpdatePositions();
    
    // Apply thermostat
    ApplyThermostat();
}

void AFEPCalculator::CalculateForces(float Lambda)
{
    // Zero forces
    for (FAtomState& Atom : CurrentState.Atoms)
    {
        Atom.Force = FVector::ZeroVector;
    }
    
    // Calculate non-bonded forces
    CalculateElectrostaticForces(Lambda);
    CalculateVanDerWaalsForces(Lambda);
    
    // Apply restraints to protein atoms
    if (Parameters.bRestrainProtein)
    {
        ApplyRestraints();
    }
}

void AFEPCalculator::UpdatePositions()
{
    float dt = Parameters.TimeStep;
    
    for (FAtomState& Atom : CurrentState.Atoms)
    {
        // Update velocity: v(t+dt/2) = v(t) + F/m * dt/2
        FVector Acceleration = Atom.Force / Atom.Mass;
        Atom.Velocity += Acceleration * (dt * 0.5f);
        
        // Update position: r(t+dt) = r(t) + v(t+dt/2) * dt
        Atom.Position += Atom.Velocity * dt;
    }
}

void AFEPCalculator::ApplyThermostat()
{
    // Berendsen thermostat for temperature control
    
    // Calculate current kinetic energy
    float KE = 0.0f;
    for (const FAtomState& Atom : CurrentState.Atoms)
    {
        KE += 0.5f * Atom.Mass * Atom.Velocity.SizeSquared();
    }
    
    CurrentState.KineticEnergy = KE;
    
    // Calculate current temperature
    // KE = (3/2) * N * kB * T
    int32 DegreesOfFreedom = 3 * CurrentState.Atoms.Num();
    float CurrentTemp = (2.0f * KE) / (DegreesOfFreedom * FEPConstants::BOLTZMANN);
    CurrentState.Temperature = CurrentTemp;
    
    // Apply velocity scaling (Berendsen coupling)
    float TargetTemp = Parameters.Temperature;
    float TauT = 0.1f; // Coupling time constant (ps)
    float Lambda = FMath::Sqrt(1.0f + (Parameters.TimeStep / TauT) * ((TargetTemp / CurrentTemp) - 1.0f));
    
    for (FAtomState& Atom : CurrentState.Atoms)
    {
        Atom.Velocity *= Lambda;
    }
}

float AFEPCalculator::CalculateTotalEnergy(float Lambda)
{
    float Elec = CalculateElectrostaticEnergy(Lambda);
    float VdW = CalculateVanDerWaalsEnergy(Lambda);
    float Solv = Parameters.bCalculateSolvation ? CalculateSolvationEnergy(Lambda) : 0.0f;
    
    float PE = Elec + VdW + Solv;
    CurrentState.PotentialEnergy = PE;
    
    return PE + CurrentState.KineticEnergy;
}

float AFEPCalculator::CalculateElectrostaticEnergy(float Lambda)
{
    float Energy = 0.0f;
    int32 NumAtoms = CurrentState.Atoms.Num();
    
    for (int32 i = 0; i < NumAtoms; ++i)
    {
        const FAtomState& Atom1 = CurrentState.Atoms[i];
        
        for (int32 j = i + 1; j < NumAtoms; ++j)
        {
            const FAtomState& Atom2 = CurrentState.Atoms[j];
            
            float r = FVector::Dist(Atom1.Position, Atom2.Position);
            
            if (r < Parameters.CutoffDistance)
            {
                // Scale ligand interactions by lambda
                float ScaleFactor = 1.0f;
                if (Atom1.bIsLigandAtom || Atom2.bIsLigandAtom)
                {
                    if (Parameters.bUseSoftCore)
                    {
                        Energy += SoftCoreCoulomb(r, Atom1.Charge, Atom2.Charge, Lambda);
                        continue;
                    }
                    else
                    {
                        ScaleFactor = Lambda;
                    }
                }
                
                // Coulomb's law: E = k * q1 * q2 / (ε * r)
                float E = (FEPConstants::COULOMB * Atom1.Charge * Atom2.Charge) / 
                          (Parameters.DielectricConstant * r);
                Energy += ScaleFactor * E;
            }
        }
    }
    
    if (Parameters.bCalculateComponents)
    {
        LastResult.ElectrostaticContribution = Energy;
    }
    
    return Energy;
}

float AFEPCalculator::CalculateVanDerWaalsEnergy(float Lambda)
{
    float Energy = 0.0f;
    int32 NumAtoms = CurrentState.Atoms.Num();
    
    for (int32 i = 0; i < NumAtoms; ++i)
    {
        const FAtomState& Atom1 = CurrentState.Atoms[i];
        
        for (int32 j = i + 1; j < NumAtoms; ++j)
        {
            const FAtomState& Atom2 = CurrentState.Atoms[j];
            
            float r = FVector::Dist(Atom1.Position, Atom2.Position);
            
            if (r < Parameters.CutoffDistance)
            {
                // Lennard-Jones parameters
                float Sigma = (Atom1.VdWRadius + Atom2.VdWRadius) * 0.5f;
                float Epsilon = FMath::Sqrt(Atom1.VdWEpsilon * Atom2.VdWEpsilon);
                
                // Scale ligand interactions by lambda
                float ScaleFactor = 1.0f;
                if (Atom1.bIsLigandAtom || Atom2.bIsLigandAtom)
                {
                    if (Parameters.bUseSoftCore)
                    {
                        Energy += SoftCoreLJ(r, Sigma, Epsilon, Lambda);
                        continue;
                    }
                    else
                    {
                        ScaleFactor = Lambda;
                    }
                }
                
                // Lennard-Jones 12-6 potential
                float Ratio = Sigma / r;
                float R6 = FMath::Pow(Ratio, 6.0f);
                float R12 = R6 * R6;
                float E = 4.0f * Epsilon * (R12 - R6);
                
                Energy += ScaleFactor * E;
            }
        }
    }
    
    if (Parameters.bCalculateComponents)
    {
        LastResult.VdWContribution = Energy;
    }
    
    return Energy;
}

float AFEPCalculator::CalculateSolvationEnergy(float Lambda)
{
    // Simplified generalized Born solvation model
    float Energy = 0.0f;
    
    // For each ligand atom, estimate solvation penalty
    for (const FAtomState& Atom : CurrentState.Atoms)
    {
        if (!Atom.bIsLigandAtom) continue;
        
        // Estimate Born radius (simplified)
        float BornRadius = Atom.VdWRadius * 1.2f;
        
        // Self-energy term
        float SelfEnergy = -(FEPConstants::COULOMB * Atom.Charge * Atom.Charge) / 
                           (2.0f * BornRadius) * (1.0f / Parameters.DielectricConstant - 1.0f / 78.5f);
        
        Energy += Lambda * SelfEnergy;
    }
    
    if (Parameters.bCalculateComponents)
    {
        LastResult.SolvationContribution = Energy;
    }
    
    return Energy;
}

float AFEPCalculator::CalculatedHdLambda(float Lambda)
{
    // Numerical derivative: dH/dλ ≈ [H(λ+δ) - H(λ-δ)] / (2δ)
    float Delta = 0.001f;
    
    float EPlus = CalculateElectrostaticEnergy(Lambda + Delta) + 
                  CalculateVanDerWaalsEnergy(Lambda + Delta);
    float EMinus = CalculateElectrostaticEnergy(Lambda - Delta) + 
                   CalculateVanDerWaalsEnergy(Lambda - Delta);
    
    if (Parameters.bCalculateSolvation)
    {
        EPlus += CalculateSolvationEnergy(Lambda + Delta);
        EMinus += CalculateSolvationEnergy(Lambda - Delta);
    }
    
    return (EPlus - EMinus) / (2.0f * Delta);
}

void AFEPCalculator::CalculateElectrostaticForces(float Lambda)
{
    int32 NumAtoms = CurrentState.Atoms.Num();
    
    for (int32 i = 0; i < NumAtoms; ++i)
    {
        FAtomState& Atom1 = CurrentState.Atoms[i];
        
        for (int32 j = i + 1; j < NumAtoms; ++j)
        {
            FAtomState& Atom2 = CurrentState.Atoms[j];
            
            FVector Delta = Atom2.Position - Atom1.Position;
            float r = Delta.Size();
            
            if (r < Parameters.CutoffDistance && r > 0.01f)
            {
                float ScaleFactor = 1.0f;
                if (Atom1.bIsLigandAtom || Atom2.bIsLigandAtom)
                {
                    ScaleFactor = Lambda;
                }
                
                // F = -dE/dr * (r_vec / r)
                float ForceMag = ScaleFactor * (FEPConstants::COULOMB * Atom1.Charge * Atom2.Charge) / 
                                (Parameters.DielectricConstant * r * r);
                
                FVector Force = (Delta / r) * ForceMag;
                
                Atom1.Force -= Force;
                Atom2.Force += Force;
            }
        }
    }
}

void AFEPCalculator::CalculateVanDerWaalsForces(float Lambda)
{
    int32 NumAtoms = CurrentState.Atoms.Num();
    
    for (int32 i = 0; i < NumAtoms; ++i)
    {
        FAtomState& Atom1 = CurrentState.Atoms[i];
        
        for (int32 j = i + 1; j < NumAtoms; ++j)
        {
            FAtomState& Atom2 = CurrentState.Atoms[j];
            
            FVector Delta = Atom2.Position - Atom1.Position;
            float r = Delta.Size();
            
            if (r < Parameters.CutoffDistance && r > 0.01f)
            {
                float Sigma = (Atom1.VdWRadius + Atom2.VdWRadius) * 0.5f;
                float Epsilon = FMath::Sqrt(Atom1.VdWEpsilon * Atom2.VdWEpsilon);
                
                float ScaleFactor = 1.0f;
                if (Atom1.bIsLigandAtom || Atom2.bIsLigandAtom)
                {
                    ScaleFactor = Lambda;
                }
                
                // F = -dE/dr for LJ potential
                float Ratio = Sigma / r;
                float R6 = FMath::Pow(Ratio, 6.0f);
                float R12 = R6 * R6;
                
                float ForceMag = ScaleFactor * 24.0f * Epsilon * (2.0f * R12 - R6) / r;
                
                FVector Force = (Delta / r) * ForceMag;
                
                Atom1.Force -= Force;
                Atom2.Force += Force;
            }
        }
    }
}

void AFEPCalculator::ApplyRestraints()
{
    int32 LigandAtomCount = 0;
    for (const FAtomState& Atom : CurrentState.Atoms)
    {
        if (Atom.bIsLigandAtom) LigandAtomCount++;
    }
    
    int32 ProteinIndex = 0;
    for (int32 i = LigandAtomCount; i < CurrentState.Atoms.Num(); ++i)
    {
        FAtomState& Atom = CurrentState.Atoms[i];
        
        if (ProteinIndex < InitialProteinPositions.Num())
        {
            // Harmonic restraint: F = -k * (r - r0)
            FVector Displacement = Atom.Position - InitialProteinPositions[ProteinIndex];
            FVector RestraintForce = -Parameters.ProteinRestraintForce * Displacement;
            
            Atom.Force += RestraintForce;
            ProteinIndex++;
        }
    }
}

float AFEPCalculator::SoftCoreLJ(float r, float sigma, float epsilon, float lambda)
{
    float Alpha = Parameters.SoftCoreAlpha;
    float SigmaSq = sigma * sigma;
    float rSoft = FMath::Pow(Alpha * SigmaSq * (1.0f - lambda) + r * r, 0.5f);
    
    float Ratio = sigma / rSoft;
    float R6 = FMath::Pow(Ratio, 6.0f);
    float R12 = R6 * R6;
    
    return lambda * 4.0f * epsilon * (R12 - R6);
}

float AFEPCalculator::SoftCoreCoulomb(float r, float q1, float q2, float lambda)
{
    float Alpha = Parameters.SoftCoreAlpha;
    float rSoft = FMath::Pow(Alpha * (1.0f - lambda) + r * r, 0.5f);
    
    return lambda * (FEPConstants::COULOMB * q1 * q2) / (Parameters.DielectricConstant * rSoft);
}

float AFEPCalculator::IntegrateFreeEnergy(const TArray<FFEPLambdaWindow>& Windows)
{
    // Trapezoidal rule integration: ΔG = ∫(dH/dλ)dλ
    float DeltaG = 0.0f;
    float ErrorSum = 0.0f;
    
    for (int32 i = 0; i < Windows.Num() - 1; ++i)
    {
        float dLambda = Windows[i + 1].Lambda - Windows[i].Lambda;
        float AvgdHdL = (Windows[i].dHdLambda + Windows[i + 1].dHdLambda) * 0.5f;
        
        DeltaG += AvgdHdL * dLambda;
        
        // Propagate errors
        float AvgError = (Windows[i].StandardError + Windows[i + 1].StandardError) * 0.5f;
        ErrorSum += (AvgError * dLambda) * (AvgError * dLambda);
    }
    
    LastResult.StandardError = FMath::Sqrt(ErrorSum);
    
    return DeltaG;
}

float AFEPCalculator::CalculateBindingAffinityFromDeltaG(float DeltaG)
{
    // ΔG = -RT ln(K_a) = RT ln(K_d)
    // K_d = exp(ΔG / RT)
    
    float RT = FEPConstants::GAS_CONSTANT * Parameters.Temperature; // kcal/mol
    float Kd = FMath::Exp(DeltaG / RT); // Dissociation constant in M
    
    // Convert to nM
    float KdNanoMolar = Kd * 1.0e9f;
    
    return KdNanoMolar;
}

void AFEPCalculator::SetAtomParameters(FAtomState& Atom, const FString& Element)
{
    Atom.Mass = GetAtomicMass(Element);
    Atom.VdWRadius = GetVdWRadius(Element);
    Atom.VdWEpsilon = GetVdWEpsilon(Element);
    
    // Simplified charge assignment (should use force field)
    if (Element == TEXT("N") || Element == TEXT("O"))
    {
        Atom.Charge = -0.5f;
    }
    else if (Element == TEXT("H"))
    {
        Atom.Charge = 0.3f;
    }
    else if (Element == TEXT("C"))
    {
        Atom.Charge = 0.1f;
    }
    else
    {
        Atom.Charge = 0.0f;
    }
}

float AFEPCalculator::GetAtomicMass(const FString& Element)
{
    static TMap<FString, float> Masses = {
        {TEXT("H"), 1.008f},
        {TEXT("C"), 12.011f},
        {TEXT("N"), 14.007f},
        {TEXT("O"), 15.999f},
        {TEXT("F"), 18.998f},
        {TEXT("P"), 30.974f},
        {TEXT("S"), 32.065f},
        {TEXT("Cl"), 35.453f},
        {TEXT("Br"), 79.904f},
        {TEXT("I"), 126.904f}
    };
    
    return Masses.Contains(Element) ? Masses[Element] : 12.0f;
}

float AFEPCalculator::GetVdWRadius(const FString& Element)
{
    static TMap<FString, float> Radii = {
        {TEXT("H"), 1.2f},
        {TEXT("C"), 1.7f},
        {TEXT("N"), 1.55f},
        {TEXT("O"), 1.52f},
        {TEXT("F"), 1.47f},
        {TEXT("P"), 1.8f},
        {TEXT("S"), 1.8f},
        {TEXT("Cl"), 1.75f},
        {TEXT("Br"), 1.85f}
    };
    
    return Radii.Contains(Element) ? Radii[Element] : 1.7f;
}

float AFEPCalculator::GetVdWEpsilon(const FString& Element)
{
    // OPLS-AA epsilon values (kcal/mol)
    static TMap<FString, float> Epsilons = {
        {TEXT("H"), 0.03f},
        {TEXT("C"), 0.07f},
        {TEXT("N"), 0.17f},
        {TEXT("O"), 0.21f},
        {TEXT("F"), 0.061f},
        {TEXT("P"), 0.2f},
        {TEXT("S"), 0.25f},
        {TEXT("Cl"), 0.265f}
    };
    
    return Epsilons.Contains(Element) ? Epsilons[Element] : 0.1f;
}

float AFEPCalculator::CalculateBlockAverageError(const TArray<float>& Data)
{
    if (Data.Num() < 10) return 0.0f;
    
    // Block averaging for error estimation
    int32 BlockSize = FMath::Max(1, Data.Num() / 10);
    TArray<float> BlockAverages;
    
    for (int32 i = 0; i < Data.Num(); i += BlockSize)
    {
        float Sum = 0.0f;
        int32 Count = 0;
        
        for (int32 j = i; j < FMath::Min(i + BlockSize, Data.Num()); ++j)
        {
            Sum += Data[j];
            Count++;
        }
        
        if (Count > 0)
        {
            BlockAverages.Add(Sum / Count);
        }
    }
    
    // Calculate standard deviation of block averages
    float Mean = 0.0f;
    for (float Val : BlockAverages)
    {
        Mean += Val;
    }
    Mean /= BlockAverages.Num();
    
    float Variance = 0.0f;
    for (float Val : BlockAverages)
    {
        float Diff = Val - Mean;
        Variance += Diff * Diff;
    }
    Variance /= (BlockAverages.Num() - 1);
    
    // Standard error
    return FMath::Sqrt(Variance / BlockAverages.Num());
}

// Visualization methods
void AFEPCalculator::VisualizeEnergyLandscape(const FFEPResult& Result)
{
    ClearVisualization();
    
    UE_LOG(LogTemp, Log, TEXT("FEPCalculator: Visualizing energy landscape..."));
    
    // Create visualization showing lambda windows
    // This would create meshes showing the free energy profile
    // Implementation depends on your visualization preferences
}

void AFEPCalculator::ClearVisualization()
{
    for (UStaticMeshComponent* Mesh : VisualizationMeshes)
    {
        if (Mesh)
        {
            Mesh->DestroyComponent();
        }
    }
    VisualizationMeshes.Empty();
}

float AFEPCalculator::CalculateBondEnergy()
{
    // TODO: Implement bond stretching energy
    // E = k * (r - r0)^2
    return 0.0f;
}

float AFEPCalculator::CalculateAngleEnergy()
{
    // TODO: Implement angle bending energy
    // E = k * (θ - θ0)^2
    return 0.0f;
}

float AFEPCalculator::CalculateDihedralEnergy()
{
    // TODO: Implement dihedral torsion energy
    // E = sum[Vn/2 * (1 + cos(n*φ - γ))]
    return 0.0f;
}

void AFEPCalculator::CalculateBondForces()
{
    // TODO: Implement bonded force calculations
}
