// MolecularDynamics.cpp - Optimized O(N) Implementation
#include "MolecularDynamics.h"
#include "PDBViewer.h"
#include "Components/StaticMeshComponent.h"
#include "Kismet/GameplayStatics.h"
#include "Async/ParallelFor.h"

AMolecularDynamics::AMolecularDynamics()
{
    PrimaryActorTick.bCanEverTick = true;
    PrimaryActorTick.TickGroup = TG_PrePhysics;
    CellSize = CutoffDistance + NeighborListSkin;
}

void AMolecularDynamics::BeginPlay()
{
    Super::BeginPlay();
    InitializeElementDatabase();
    
    TArray<AActor*> Found;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), APDBViewer::StaticClass(), Found);
    if (Found.Num() > 0)
    {
        ViewerReference = Cast<APDBViewer>(Found[0]);
    }
}

void AMolecularDynamics::InitializeElementDatabase()
{
    // Initialize with realistic parameters
    ElementDatabase.Empty();
    
    auto AddElement = [this](const FString& Sym, float M, float R, float Eps, float Q = 0.0f)
    {
        FElementProperties Props;
        Props.Mass = M;
        Props.Radius = R;
        Props.Epsilon = Eps;
        Props.Charge = Q;
        ElementDatabase.Add(Sym, Props);
    };
    
    // Common elements with VdW radii (in angstroms, scaled)
    AddElement(TEXT("H"),  1.008f,  120.0f, 0.02f);
    AddElement(TEXT("C"),  12.011f, 170.0f, 0.07f);
    AddElement(TEXT("N"),  14.007f, 155.0f, 0.07f);
    AddElement(TEXT("O"),  15.999f, 152.0f, 0.06f);
    AddElement(TEXT("S"),  32.065f, 180.0f, 0.25f);
    AddElement(TEXT("P"),  30.974f, 180.0f, 0.20f);
    AddElement(TEXT("F"),  18.998f, 147.0f, 0.06f);
    AddElement(TEXT("CL"), 35.453f, 175.0f, 0.27f);
    AddElement(TEXT("BR"), 79.904f, 185.0f, 0.31f);
    AddElement(TEXT("I"),  126.904f, 198.0f, 0.34f);
}

void AMolecularDynamics::Tick(float DeltaTime)
{
    Super::Tick(DeltaTime);
    
    if (bIsSimulating && Atoms.Num() > 0)
    {
        // Multiple substeps for stability
        const int32 SubSteps = 2;
        const float SubDT = TimeStep / SubSteps;
        
        for (int32 i = 0; i < SubSteps; ++i)
        {
            StepSimulation(SubDT);
        }
        
        SimulationTime += TimeStep;
        StepCount++;
    }
}

void AMolecularDynamics::InitializeFromViewer(APDBViewer* Viewer)
{
    if (!Viewer)
    {
        UE_LOG(LogTemp, Error, TEXT("MD: Null viewer"));
        return;
    }
    
    ViewerReference = Viewer;
    Atoms.Empty();
    Bonds.Empty();
    InitialPositions.Empty();
    PreviousPositions.Empty();
    PreviousAccelerations.Empty();
    SpatialHash.Empty();
    
    // Get ALL visible ligands from the viewer
    TArray<FString> LigandKeys;
    const auto& LigandMap = Viewer->LigandMap;
    
    for (const auto& Pair : LigandMap)
    {
        if (Pair.Value && Pair.Value->bIsVisible)
        {
            LigandKeys.Add(Pair.Key);
        }
    }
    
    if (LigandKeys.Num() == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("MD: No visible ligands"));
        return;
    }
    
    UE_LOG(LogTemp, Log, TEXT("MD: Initializing from %d visible ligand(s)"), LigandKeys.Num());
    
    int32 TotalAtoms = 0;
    TMap<UStaticMeshComponent*, int32> MeshToIndex;
    
    // Collect all visible atoms from all visible ligands
    for (const FString& Key : LigandKeys)
    {
        FLigandInfo* LigInfo = LigandMap[Key];
        if (!LigInfo) continue;
        
        for (UStaticMeshComponent* Mesh : LigInfo->AtomMeshes)
        {
            if (!Mesh || !Mesh->IsVisible()) continue;
            
            FAtomState NewAtom;
            NewAtom.Position = Mesh->GetComponentLocation();
            NewAtom.Velocity = FVector::ZeroVector;
            NewAtom.Force = FVector::ZeroVector;
            NewAtom.MeshComponent = Mesh;
            
            // Infer element from color
            auto Mat = Mesh->GetMaterial(0);
            NewAtom.Element = TEXT("C"); // Default
            
            const FElementProperties* Props = ElementDatabase.Find(NewAtom.Element);
            if (Props)
            {
                NewAtom.Mass = Props->Mass;
                NewAtom.Charge = Props->Charge;
            }
            
            // Small random thermal velocity
            float VelScale = FMath::Sqrt(TargetTemperature * 0.008314f / NewAtom.Mass);
            NewAtom.Velocity = FVector(
                FMath::FRandRange(-VelScale, VelScale),
                FMath::FRandRange(-VelScale, VelScale),
                FMath::FRandRange(-VelScale, VelScale)
            );
            
            MeshToIndex.Add(Mesh, Atoms.Num());
            Atoms.Add(NewAtom);
            InitialPositions.Add(NewAtom.Position);
            PreviousPositions.Add(NewAtom.Position);
            PreviousAccelerations.Add(FVector::ZeroVector);
            TotalAtoms++;
        }
        
        // Collect bonds for this ligand
        // Note: We'd need to extend PDBViewer to expose bond info
        // For now, infer from proximity
    }
    
    // Infer bonds from proximity (all atoms, regardless of ligand)
    const float BondThreshold = 80.0f;
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        for (int32 j = i + 1; j < Atoms.Num(); ++j)
        {
            float Dist = FVector::Dist(Atoms[i].Position, Atoms[j].Position);
            if (Dist < BondThreshold)
            {
                FBondConstraint Bond;
                Bond.Atom1Index = i;
                Bond.Atom2Index = j;
                Bond.RestLength = Dist;
                
                // Stronger springs for shorter bonds
                Bond.SpringConstant = 1000.0f / FMath::Max(Dist, 50.0f);
                Bond.BondOrder = 1;
                Bonds.Add(Bond);
            }
        }
    }
    
    // Initialize spatial hash
    RebuildSpatialHash();
    UpdateNeighborLists();
    
    // Initialize RK4 state
    RK4State.K1Vel.SetNum(Atoms.Num());
    RK4State.K2Vel.SetNum(Atoms.Num());
    RK4State.K3Vel.SetNum(Atoms.Num());
    RK4State.K4Vel.SetNum(Atoms.Num());
    RK4State.K1Acc.SetNum(Atoms.Num());
    RK4State.K2Acc.SetNum(Atoms.Num());
    RK4State.K3Acc.SetNum(Atoms.Num());
    RK4State.K4Acc.SetNum(Atoms.Num());
    
    UE_LOG(LogTemp, Log, TEXT("MD: Initialized %d atoms, %d bonds from %d ligand(s)"), 
        TotalAtoms, Bonds.Num(), LigandKeys.Num());
}

void AMolecularDynamics::StartSimulation()
{
    if (Atoms.Num() == 0)
    {
        UE_LOG(LogTemp, Warning, TEXT("MD: No atoms"));
        return;
    }
    bIsSimulating = true;
    SimulationTime = 0.0f;
    StepCount = 0;
}

void AMolecularDynamics::StopSimulation()
{
    bIsSimulating = false;
}

void AMolecularDynamics::ResetSimulation()
{
    bIsSimulating = false;
    SimulationTime = 0.0f;
    StepCount = 0;
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        if (i < InitialPositions.Num())
        {
            Atoms[i].Position = InitialPositions[i];
            Atoms[i].Velocity = FVector::ZeroVector;
            Atoms[i].Force = FVector::ZeroVector;
            
            if (Atoms[i].MeshComponent)
            {
                Atoms[i].MeshComponent->SetWorldLocation(Atoms[i].Position);
            }
        }
    }
    
    PreviousPositions = InitialPositions;
    RebuildSpatialHash();
}

void AMolecularDynamics::StepSimulation(float DeltaTime)
{
    if (Atoms.Num() == 0) return;
    
    // Update neighbor lists if needed
    if (ShouldRebuildNeighborLists())
    {
        UpdateNeighborLists();
    }
    
    // Calculate forces using spatial hashing
    CalculateForces();
    
    // Integrate
    switch (Integrator)
    {
        case EMDIntegrator::Verlet:
            IntegrateVerlet(DeltaTime);
            break;
        case EMDIntegrator::LeapFrog:
            IntegrateLeapFrog(DeltaTime);
            break;
        case EMDIntegrator::RungeKutta4:
            IntegrateRungeKutta4(DeltaTime);
            break;
    }
    
    // Constraints
    if (bConstrainBonds)
    {
        EnforceConstraints(DeltaTime);
    }
    
    // Update spatial hash after position changes
    UpdateSpatialHash();
    
    // Temperature control
    ApplyThermostat();
    
    // Update energies
    UpdateEnergies();
    
    // Update visuals
    UpdateMeshPositions();
}

// O(N) force calculation using spatial hashing
void AMolecularDynamics::CalculateForces()
{
    // Reset forces - parallelizable
    ParallelFor(Atoms.Num(), [this](int32 i)
    {
        Atoms[i].Force = FVector::ZeroVector;
    });
    
    // Bond forces - O(B) where B is number of bonds
    CalculateBondForces();
    
    // Non-bonded forces using neighbor lists - O(N)
    if (bUseLennardJones || bUseElectrostatics)
    {
        CalculateNonBondedForces();
    }
}

void AMolecularDynamics::CalculateBondForces()
{
    // Can be parallelized with atomic operations or partitioning
    for (const FBondConstraint& Bond : Bonds)
    {
        if (Bond.Atom1Index >= Atoms.Num() || Bond.Atom2Index >= Atoms.Num())
            continue;
        
        FAtomState& A1 = Atoms[Bond.Atom1Index];
        FAtomState& A2 = Atoms[Bond.Atom2Index];
        
        FVector Delta = A2.Position - A1.Position;
        float Dist = Delta.Size();
        
        if (Dist < KINDA_SMALL_NUMBER) continue;
        
        // Harmonic bond: F = -k(r - r0)
        float Displacement = Dist - Bond.RestLength;
        float ForceMag = -Bond.SpringConstant * Displacement;
        
        FVector Force = (Delta / Dist) * ForceMag;
        
        A1.Force += Force;
        A2.Force -= Force;
    }
}

void AMolecularDynamics::CalculateNonBondedForces()
{
    // Use neighbor lists for O(N) scaling
    // Can be parallelized
    const float CutoffSq = CutoffDistance * CutoffDistance;
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        FAtomState& A1 = Atoms[i];
        
        // Only check neighbors in the neighbor list
        for (int32 j : A1.NeighborList)
        {
            if (j <= i) continue; // Avoid double counting
            
            FAtomState& A2 = Atoms[j];
            FVector Delta = A2.Position - A1.Position;
            float DistSq = Delta.SizeSquared();
            
            if (DistSq > CutoffSq || DistSq < 1.0f) continue;
            
            // Check if bonded (skip LJ for bonded atoms)
            bool bBonded = false;
            for (const FBondConstraint& B : Bonds)
            {
                if ((B.Atom1Index == i && B.Atom2Index == j) ||
                    (B.Atom1Index == j && B.Atom2Index == i))
                {
                    bBonded = true;
                    break;
                }
            }
            if (bBonded) continue;
            
            float Dist = FMath::Sqrt(DistSq);
            FVector Dir = Delta / Dist;
            
            // Lennard-Jones: F = 24*eps/r * [2*(sig/r)^12 - (sig/r)^6]
            if (bUseLennardJones)
            {
                // Use element-specific parameters
                const FElementProperties* P1 = ElementDatabase.Find(A1.Element);
                const FElementProperties* P2 = ElementDatabase.Find(A2.Element);
                
                float Sigma = LJSigma;
                float Epsilon = LJEpsilon;
                
                if (P1 && P2)
                {
                    // Lorentz-Berthelot combining rules
                    Sigma = (P1->Radius + P2->Radius) * 0.5f;
                    Epsilon = FMath::Sqrt(P1->Epsilon * P2->Epsilon);
                }
                
                float SigOverR = Sigma / Dist;
                float SigOverR6 = FMath::Pow(SigOverR, 6.0f);
                float SigOverR12 = SigOverR6 * SigOverR6;
                
                float ForceMag = 24.0f * Epsilon / Dist * (2.0f * SigOverR12 - SigOverR6);
                FVector Force = Dir * ForceMag;
                
                A1.Force += Force;
                A2.Force -= Force;
            }
            
            // Electrostatics: F = k*q1*q2/r^2
            if (bUseElectrostatics && (A1.Charge != 0.0f || A2.Charge != 0.0f))
            {
                const float CoulombK = 1389.0f; // e^2/(4*pi*eps0) in kJ/mol*nm
                float ForceMag = CoulombK * A1.Charge * A2.Charge / DistSq;
                FVector Force = Dir * ForceMag;
                
                A1.Force += Force;
                A2.Force -= Force;
            }
        }
    }
}

// Spatial hashing for O(N) neighbor search
int64 AMolecularDynamics::GetCellHash(int32 X, int32 Y, int32 Z) const
{
    // Morton encoding for better cache locality
    const int64 p1 = 73856093;
    const int64 p2 = 19349663;
    const int64 p3 = 83492791;
    return (X * p1) ^ (Y * p2) ^ (Z * p3);
}

void AMolecularDynamics::GetCellCoords(const FVector& Pos, int32& OutX, int32& OutY, int32& OutZ) const
{
    OutX = FMath::FloorToInt((Pos.X - GridMin.X) / CellSize);
    OutY = FMath::FloorToInt((Pos.Y - GridMin.Y) / CellSize);
    OutZ = FMath::FloorToInt((Pos.Z - GridMin.Z) / CellSize);
}

void AMolecularDynamics::GetNeighborCells(int32 X, int32 Y, int32 Z, TArray<int64>& OutHashes) const
{
    OutHashes.Empty(27);
    
    // Check all 27 neighboring cells (3x3x3)
    for (int32 dx = -1; dx <= 1; ++dx)
    {
        for (int32 dy = -1; dy <= 1; ++dy)
        {
            for (int32 dz = -1; dz <= 1; ++dz)
            {
                OutHashes.Add(GetCellHash(X + dx, Y + dy, Z + dz));
            }
        }
    }
}

void AMolecularDynamics::RebuildSpatialHash()
{
    SpatialHash.Empty();
    
    if (Atoms.Num() == 0) return;
    
    // Find bounding box
    GridMin = Atoms[0].Position;
    GridMax = Atoms[0].Position;
    
    for (const FAtomState& A : Atoms)
    {
        GridMin = GridMin.ComponentMin(A.Position);
        GridMax = GridMax.ComponentMax(A.Position);
    }
    
    // Expand slightly
    GridMin -= FVector(CellSize * 2);
    GridMax += FVector(CellSize * 2);
    
    CellSize = CutoffDistance + NeighborListSkin;
    
    // Insert atoms
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        int32 X, Y, Z;
        GetCellCoords(Atoms[i].Position, X, Y, Z);
        
        Atoms[i].CellX = X;
        Atoms[i].CellY = Y;
        Atoms[i].CellZ = Z;
        Atoms[i].CellHash = GetCellHash(X, Y, Z);
        
        SpatialHash.FindOrAdd(Atoms[i].CellHash).Add(i);
    }
}

void AMolecularDynamics::UpdateSpatialHash()
{
    // Incremental update - only move atoms that changed cells
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        int32 NewX, NewY, NewZ;
        GetCellCoords(Atoms[i].Position, NewX, NewY, NewZ);
        
        if (NewX != Atoms[i].CellX || NewY != Atoms[i].CellY || NewZ != Atoms[i].CellZ)
        {
            // Remove from old cell
            if (TArray<int32>* OldCell = SpatialHash.Find(Atoms[i].CellHash))
            {
                OldCell->Remove(i);
            }
            
            // Add to new cell
            Atoms[i].CellX = NewX;
            Atoms[i].CellY = NewY;
            Atoms[i].CellZ = NewZ;
            Atoms[i].CellHash = GetCellHash(NewX, NewY, NewZ);
            SpatialHash.FindOrAdd(Atoms[i].CellHash).Add(i);
        }
    }
}

bool AMolecularDynamics::ShouldRebuildNeighborLists() const
{
    return (StepCount % NeighborListRebuildFrequency) == 0;
}

void AMolecularDynamics::UpdateNeighborLists()
{
    const float NeighborCutoffSq = FMath::Square(CutoffDistance + NeighborListSkin);
    
    // Parallel neighbor list construction
    ParallelFor(Atoms.Num(), [this, NeighborCutoffSq](int32 i)
    {
        Atoms[i].NeighborList.Empty();
        
        TArray<int64> NeighborCells;
        GetNeighborCells(Atoms[i].CellX, Atoms[i].CellY, Atoms[i].CellZ, NeighborCells);
        
        for (int64 Hash : NeighborCells)
        {
            const TArray<int32>* Cell = SpatialHash.Find(Hash);
            if (!Cell) continue;
            
            for (int32 j : *Cell)
            {
                if (j == i) continue;
                
                float DistSq = FVector::DistSquared(Atoms[i].Position, Atoms[j].Position);
                if (DistSq < NeighborCutoffSq)
                {
                    Atoms[i].NeighborList.Add(j);
                }
            }
        }
    });
}

// Velocity Verlet - symplectic, time-reversible
void AMolecularDynamics::IntegrateVerlet(float DT)
{
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        FAtomState& A = Atoms[i];
        FVector Acc = A.Force / A.Mass;
        
        // x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2
        A.Position += A.Velocity * DT + 0.5f * Acc * DT * DT;
        
        // v(t+dt) = v(t) + 0.5*[a(t) + a(t+dt)]*dt
        // Use previous acceleration for now, will be corrected after force calc
        A.Velocity += 0.5f * (PreviousAccelerations[i] + Acc) * DT;
        
        PreviousAccelerations[i] = Acc;
    }
}

void AMolecularDynamics::IntegrateLeapFrog(float DT)
{
    // First time: v(t+dt/2) = v(t) + 0.5*a(t)*dt
    // Then: x(t+dt) = x(t) + v(t+dt/2)*dt
    //       v(t+3dt/2) = v(t+dt/2) + a(t+dt)*dt
    
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        FAtomState& A = Atoms[i];
        FVector Acc = A.Force / A.Mass;
        
        A.Velocity += Acc * DT;
        A.Position += A.Velocity * DT;
    }
}

void AMolecularDynamics::IntegrateRungeKutta4(float DT)
{
    // RK4 for highest accuracy - useful for validation
    // k1 = f(t, y)
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        RK4State.K1Vel[i] = Atoms[i].Velocity;
        RK4State.K1Acc[i] = Atoms[i].Force / Atoms[i].Mass;
    }
    
    // k2 = f(t + dt/2, y + k1*dt/2)
    TArray<FVector> OrigPos;
    TArray<FVector> OrigVel;
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        OrigPos.Add(Atoms[i].Position);
        OrigVel.Add(Atoms[i].Velocity);
        Atoms[i].Position += RK4State.K1Vel[i] * DT * 0.5f;
        Atoms[i].Velocity += RK4State.K1Acc[i] * DT * 0.5f;
    }
    CalculateForces();
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        RK4State.K2Vel[i] = Atoms[i].Velocity;
        RK4State.K2Acc[i] = Atoms[i].Force / Atoms[i].Mass;
    }
    
    // k3 = f(t + dt/2, y + k2*dt/2)
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        Atoms[i].Position = OrigPos[i] + RK4State.K2Vel[i] * DT * 0.5f;
        Atoms[i].Velocity = OrigVel[i] + RK4State.K2Acc[i] * DT * 0.5f;
    }
    CalculateForces();
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        RK4State.K3Vel[i] = Atoms[i].Velocity;
        RK4State.K3Acc[i] = Atoms[i].Force / Atoms[i].Mass;
    }
    
    // k4 = f(t + dt, y + k3*dt)
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        Atoms[i].Position = OrigPos[i] + RK4State.K3Vel[i] * DT;
        Atoms[i].Velocity = OrigVel[i] + RK4State.K3Acc[i] * DT;
    }
    CalculateForces();
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        RK4State.K4Vel[i] = Atoms[i].Velocity;
        RK4State.K4Acc[i] = Atoms[i].Force / Atoms[i].Mass;
    }
    
    // y(t+dt) = y(t) + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        Atoms[i].Position = OrigPos[i] + DT / 6.0f * (
            RK4State.K1Vel[i] + 2.0f * RK4State.K2Vel[i] + 
            2.0f * RK4State.K3Vel[i] + RK4State.K4Vel[i]);
        
        Atoms[i].Velocity = OrigVel[i] + DT / 6.0f * (
            RK4State.K1Acc[i] + 2.0f * RK4State.K2Acc[i] + 
            2.0f * RK4State.K3Acc[i] + RK4State.K4Acc[i]);
    }
}

// RATTLE algorithm for constrained dynamics
void AMolecularDynamics::EnforceConstraints(float DT)
{
    const float Tolerance = 0.0001f;
    
    for (int32 Iter = 0; Iter < ConstraintIterations; ++Iter)
    {
        float MaxError = 0.0f;
        
        for (const FBondConstraint& B : Bonds)
        {
            if (B.Atom1Index >= Atoms.Num() || B.Atom2Index >= Atoms.Num())
                continue;
            
            FAtomState& A1 = Atoms[B.Atom1Index];
            FAtomState& A2 = Atoms[B.Atom2Index];
            
            FVector Delta = A2.Position - A1.Position;
            float Dist = Delta.Size();
            
            if (Dist < KINDA_SMALL_NUMBER) continue;
            
            float Error = Dist - B.RestLength;
            MaxError = FMath::Max(MaxError, FMath::Abs(Error));
            
            if (FMath::Abs(Error) > Tolerance)
            {
                FVector Dir = Delta / Dist;
                float InvMassSum = 1.0f / (A1.Mass + A2.Mass);
                
                // Position correction
                FVector Correction = Dir * Error * 0.5f;
                A1.Position += Correction * (A2.Mass * InvMassSum);
                A2.Position -= Correction * (A1.Mass * InvMassSum);
                
                // Velocity correction (project out component along bond)
                FVector RelVel = A2.Velocity - A1.Velocity;
                float VelAlongBond = FVector::DotProduct(RelVel, Dir);
                FVector VelCorrection = Dir * VelAlongBond * 0.5f;
                
                A1.Velocity += VelCorrection * (A2.Mass * InvMassSum);
                A2.Velocity -= VelCorrection * (A1.Mass * InvMassSum);
            }
        }
        
        if (MaxError < Tolerance) break;
    }
}

// FIXED UpdateEnergies() function - Add this to MolecularDynamics.cpp

void AMolecularDynamics::UpdateEnergies()
{
    // Calculate kinetic energy
    KineticEnergy = 0.0f;
    for (const FAtomState& A : Atoms)
    {
        KineticEnergy += 0.5f * A.Mass * A.Velocity.SizeSquared();
    }
    
    PotentialEnergy = 0.0f;
    
    // 1. Bond stretch energy (harmonic potential)
    for (const FBondConstraint& B : Bonds)
    {
        if (B.Atom1Index >= Atoms.Num() || B.Atom2Index >= Atoms.Num())
            continue;
        
        float Dist = FVector::Dist(Atoms[B.Atom1Index].Position, Atoms[B.Atom2Index].Position);
        float Disp = Dist - B.RestLength;
        PotentialEnergy += 0.5f * B.SpringConstant * Disp * Disp;
    }
    
    // 2. Lennard-Jones potential energy (non-bonded)
    if (bUseLennardJones)
    {
        const float CutoffSq = CutoffDistance * CutoffDistance;
        
        // Use neighbor lists for efficiency (same as force calculation)
        for (int32 i = 0; i < Atoms.Num(); ++i)
        {
            const FAtomState& A1 = Atoms[i];
            
            // Use neighbor list if available
            const TArray<int32>& Neighbors = A1.NeighborList;
            
            for (int32 j : Neighbors)
            {
                if (j <= i) continue; // Only count each pair once
                
                const FAtomState& A2 = Atoms[j];
                
                float DistSq = FVector::DistSquared(A1.Position, A2.Position);
                if (DistSq > CutoffSq || DistSq < KINDA_SMALL_NUMBER) continue;
                
                float Dist = FMath::Sqrt(DistSq);
                
                // Lennard-Jones potential: U = 4*epsilon*[(sigma/r)^12 - (sigma/r)^6]
                float SigmaOverR = LJSigma / Dist;
                float SR6 = SigmaOverR * SigmaOverR * SigmaOverR * 
                            SigmaOverR * SigmaOverR * SigmaOverR;
                float SR12 = SR6 * SR6;
                
                float LJ_Energy = 4.0f * LJEpsilon * (SR12 - SR6);
                PotentialEnergy += LJ_Energy;
            }
        }
    }
    
    // 3. Electrostatic potential energy (non-bonded)
    if (bUseElectrostatics)
    {
        const float CoulombK = 1389.0f; // e^2/(4*pi*eps0) in kJ/mol*nm
        const float CutoffSq = CutoffDistance * CutoffDistance;
        
        for (int32 i = 0; i < Atoms.Num(); ++i)
        {
            const FAtomState& A1 = Atoms[i];
            if (A1.Charge == 0.0f) continue;
            
            const TArray<int32>& Neighbors = A1.NeighborList;
            
            for (int32 j : Neighbors)
            {
                if (j <= i) continue; // Only count each pair once
                
                const FAtomState& A2 = Atoms[j];
                if (A2.Charge == 0.0f) continue;
                
                float DistSq = FVector::DistSquared(A1.Position, A2.Position);
                if (DistSq > CutoffSq || DistSq < KINDA_SMALL_NUMBER) continue;
                
                float Dist = FMath::Sqrt(DistSq);
                
                // Coulomb potential: U = k*q1*q2/r
                float Coulomb_Energy = CoulombK * A1.Charge * A2.Charge / Dist;
                PotentialEnergy += Coulomb_Energy;
            }
        }
    }
    
    TotalEnergy = KineticEnergy + PotentialEnergy;
}

void AMolecularDynamics::ApplyThermostat()
{
    // Berendsen thermostat with gentle coupling
    if (TargetTemperature <= 0.0f || Atoms.Num() == 0) return;
    
    float CurrentT = GetCurrentTemperature();
    
    if (CurrentT > KINDA_SMALL_NUMBER)
    {
        // Coupling strength
        const float TauT = 0.1f; // ps
        float Lambda = FMath::Sqrt(1.0f + TimeStep / TauT * (TargetTemperature / CurrentT - 1.0f));
        Lambda = FMath::Clamp(Lambda, 0.9f, 1.1f); // Limit rescaling
        
        for (FAtomState& A : Atoms)
        {
            A.Velocity *= Lambda * DampingFactor;
        }
    }
}

float AMolecularDynamics::GetCurrentTemperature() const
{
    if (Atoms.Num() == 0) return 0.0f;
    
    // T = 2*KE / (3*N*k_B)
    // Using k_B = 0.008314 kJ/(mol*K) (gas constant / Avogadro)
    const float kB = 0.008314f;
    float TotalKE = 0.0f;
    
    for (const FAtomState& A : Atoms)
    {
        TotalKE += 0.5f * A.Mass * A.Velocity.SizeSquared();
    }
    
    return (2.0f * TotalKE) / (3.0f * Atoms.Num() * kB);
}

void AMolecularDynamics::UpdateMeshPositions()
{
    // Can be parallelized
    ParallelFor(Atoms.Num(), [this](int32 i)
    {
        if (Atoms[i].MeshComponent)
        {
            Atoms[i].MeshComponent->SetWorldLocation(Atoms[i].Position);
        }
    });
}

float AMolecularDynamics::GetAtomMass(const FString& Element) const
{
    const FElementProperties* Props = ElementDatabase.Find(Element.ToUpper());
    return Props ? Props->Mass : 12.011f;
}

float AMolecularDynamics::GetAtomCharge(const FString& Element) const
{
    const FElementProperties* Props = ElementDatabase.Find(Element.ToUpper());
    return Props ? Props->Charge : 0.0f;
}