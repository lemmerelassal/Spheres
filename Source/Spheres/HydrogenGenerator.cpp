// HydrogenGenerator.cpp - Implementation of hydrogen generation

#include "HydrogenGenerator.h"

TArray<TPair<FVector, int32>> FHydrogenGenerator::GenerateHydrogens(
    const TArray<FVector>& AtomPositions,
    const TArray<FString>& AtomElements,
    const TArray<TPair<int32, int32>>& BondPairs,
    const TArray<int32>& BondOrders)
{
    TArray<TPair<FVector, int32>> Hydrogens;
    
    if (AtomPositions.Num() != AtomElements.Num())
        return Hydrogens;
    
    // Build bond connectivity with bond orders
    TArray<FAtomBondInfo> Atoms;
    Atoms.SetNum(AtomPositions.Num());
    
    for (int32 i = 0; i < AtomPositions.Num(); ++i)
    {
        Atoms[i].Element = AtomElements[i];
        Atoms[i].Position = AtomPositions[i];
    }
    
    // Record bonds with their orders
    for (int32 i = 0; i < BondPairs.Num(); ++i)
    {
        const auto& Bond = BondPairs[i];
        int32 Idx1 = Bond.Key;
        int32 Idx2 = Bond.Value;
        int32 Order = BondOrders.IsValidIndex(i) ? BondOrders[i] : 1;
        
        if (Atoms.IsValidIndex(Idx1) && Atoms.IsValidIndex(Idx2))
        {
            Atoms[Idx1].BondedAtomIndices.Add(Idx2);
            Atoms[Idx1].BondOrders.Add(Order);
            Atoms[Idx1].BondCount++;
            Atoms[Idx1].TotalBondOrder += Order;
            
            Atoms[Idx2].BondedAtomIndices.Add(Idx1);
            Atoms[Idx2].BondOrders.Add(Order);
            Atoms[Idx2].BondCount++;
            Atoms[Idx2].TotalBondOrder += Order;
        }
    }
    
    // Generate hydrogens for each heavy atom
    for (int32 i = 0; i < Atoms.Num(); ++i)
    {
        const FAtomBondInfo& Atom = Atoms[i];
        
        // Skip if already hydrogen
        if (Atom.Element.Equals(TEXT("H"), ESearchCase::IgnoreCase) ||
            Atom.Element.Equals(TEXT("D"), ESearchCase::IgnoreCase))
            continue;
        
        int32 NumHydrogens = GetIdealHydrogenCount(Atom.Element, Atom.TotalBondOrder);
        if (NumHydrogens <= 0)
            continue;
        
        // Get positions of bonded atoms
        TArray<FVector> BondedPositions;
        for (int32 BondedIdx : Atom.BondedAtomIndices)
        {
            if (Atoms.IsValidIndex(BondedIdx))
                BondedPositions.Add(Atoms[BondedIdx].Position);
        }
        
        // Generate hydrogen positions
        TArray<FVector> HPositions = GenerateHydrogensForAtom(
            Atom.Position,
            Atom.Element,
            BondedPositions,
            NumHydrogens,
            Atom.TotalBondOrder);
        
        // Add to result array with parent atom index
        for (const FVector& HPos : HPositions)
        {
            Hydrogens.Add(TPair<FVector, int32>(HPos, i));
        }
    }
    
    return Hydrogens;
}

int32 FHydrogenGenerator::GetIdealHydrogenCount(const FString& Element, int32 TotalBondOrder)
{
    FString Elem = Element.ToUpper();
    
    // Calculate based on total bond order (sum of all bond orders to this atom)
    // Valence - TotalBondOrder = number of H needed
    
    if (Elem == TEXT("C"))
        return FMath::Max(0, 4 - TotalBondOrder);  // sp³: 4, sp²: 3, sp: 2
    else if (Elem == TEXT("N"))
        return FMath::Max(0, 3 - TotalBondOrder);  // Can have 0-3 H depending on bonds
    else if (Elem == TEXT("O"))
        return FMath::Max(0, 2 - TotalBondOrder);  // Can have 0-2 H
    else if (Elem == TEXT("S"))
        return FMath::Max(0, 2 - TotalBondOrder);  // Similar to oxygen
    else if (Elem == TEXT("P"))
        return FMath::Max(0, 3 - TotalBondOrder);  // Can have 0-3 H
    else if (Elem == TEXT("F") || Elem == TEXT("CL") || 
             Elem == TEXT("BR") || Elem == TEXT("I"))
        return FMath::Max(0, 1 - TotalBondOrder);  // Halogens
    
    return 0;
}

float FHydrogenGenerator::GetIdealBondLength(const FString& Element)
{
    FString Elem = Element.ToUpper();
    
    // Bond lengths in Angstroms, scaled by PDB::SCALE (50.0)
    if (Elem == TEXT("C"))
        return 1.09f * 50.0f; // C-H
    else if (Elem == TEXT("N"))
        return 1.01f * 50.0f; // N-H
    else if (Elem == TEXT("O"))
        return 0.96f * 50.0f; // O-H
    else if (Elem == TEXT("S"))
        return 1.34f * 50.0f; // S-H
    else if (Elem == TEXT("P"))
        return 1.42f * 50.0f; // P-H
    
    return 1.09f * 50.0f; // Default
}

TArray<FVector> FHydrogenGenerator::GenerateHydrogensForAtom(
    const FVector& AtomPos,
    const FString& Element,
    const TArray<FVector>& BondedPositions,
    int32 NumHydrogens,
    int32 TotalBondOrder)
{
    float BondLength = GetIdealBondLength(Element);
    int32 TotalBonds = BondedPositions.Num() + NumHydrogens;
    
    // Determine geometry based on total bond order (hybridization)
    // sp³ (tetrahedral): total bond order = 4
    // sp² (trigonal): total bond order = 3 or has double bond
    // sp (linear): total bond order = 2 or has triple bond
    
    int32 ExpectedTotalOrder = TotalBondOrder + NumHydrogens;
    
    if (ExpectedTotalOrder == 4 || (TotalBonds == 4 && TotalBondOrder <= 4))
        return GenerateTetrahedralHydrogens(AtomPos, BondedPositions, NumHydrogens, BondLength);
    else if (ExpectedTotalOrder == 3 || (TotalBonds == 3 && TotalBondOrder <= 3))
        return GenerateTrigonalHydrogens(AtomPos, BondedPositions, NumHydrogens, BondLength);
    else if (ExpectedTotalOrder == 2 || TotalBonds == 2)
        return GenerateLinearHydrogens(AtomPos, BondedPositions, NumHydrogens, BondLength);
    else if (TotalBonds == 1)
    {
        // Single hydrogen - place in arbitrary direction
        TArray<FVector> Result;
        FVector Dir = GetPerpendicularVector(FVector::UpVector);
        Result.Add(AtomPos + Dir * BondLength);
        return Result;
    }
    
    return TArray<FVector>();
}

TArray<FVector> FHydrogenGenerator::GenerateTetrahedralHydrogens(
    const FVector& Center,
    const TArray<FVector>& Existing,
    int32 NumH,
    float BondLength)
{
    TArray<FVector> Hydrogens;
    
    // Tetrahedral angles: 109.5 degrees
    const float TetAngle = FMath::DegreesToRadians(109.471f);
    
    if (Existing.Num() == 3 && NumH == 1)
    {
        // 3 bonds exist, add 1 hydrogen
        FVector V1 = (Existing[0] - Center).GetSafeNormal();
        FVector V2 = (Existing[1] - Center).GetSafeNormal();
        FVector V3 = (Existing[2] - Center).GetSafeNormal();
        
        // Tetrahedral H points opposite to centroid of 3 bonds
        FVector Dir = -(V1 + V2 + V3).GetSafeNormal();
        Hydrogens.Add(Center + Dir * BondLength);
    }
    else if (Existing.Num() == 2 && NumH == 2)
    {
        // 2 bonds exist, add 2 hydrogens
        FVector V1 = (Existing[0] - Center).GetSafeNormal();
        FVector V2 = (Existing[1] - Center).GetSafeNormal();
        
        // Get perpendicular to bond plane
        FVector Normal = FVector::CrossProduct(V1, V2).GetSafeNormal();
        if (Normal.SizeSquared() < KINDA_SMALL_NUMBER)
            Normal = GetPerpendicularVector(V1);
        
        // Bisector of existing bonds
        FVector Bisector = (V1 + V2).GetSafeNormal();
        
        // Two H's symmetric about the bisector
        FVector Perp = FVector::CrossProduct(Bisector, Normal).GetSafeNormal();
        float RotAngle = FMath::Acos(FVector::DotProduct(V1, Bisector));
        float HalfAngle = (TetAngle - RotAngle) * 0.5f;
        
        FVector Dir = -Bisector;
        FQuat Q1 = FQuat(Perp, HalfAngle);
        FQuat Q2 = FQuat(Perp, -HalfAngle);
        
        Hydrogens.Add(Center + Q1.RotateVector(Dir) * BondLength);
        Hydrogens.Add(Center + Q2.RotateVector(Dir) * BondLength);
    }
    else if (Existing.Num() == 1 && NumH == 3)
    {
        // 1 bond exists, add 3 hydrogens
        FVector V = (Existing[0] - Center).GetSafeNormal();
        FVector Perp = GetPerpendicularVector(V);
        
        // Rotate 120 degrees around the existing bond
        for (int32 i = 0; i < 3; ++i)
        {
            float Angle = (i * 120.0f) * PI / 180.0f;
            FQuat Rotation = FQuat(V, Angle);
            FVector Rotated = Rotation.RotateVector(Perp);
            
            // Tilt down by tetrahedral angle
            FVector Cross = FVector::CrossProduct(V, Rotated).GetSafeNormal();
            FQuat Tilt = FQuat(Cross, PI - TetAngle);
            FVector Dir = Tilt.RotateVector(V);
            
            Hydrogens.Add(Center + Dir * BondLength);
        }
    }
    else if (Existing.Num() == 0 && NumH == 4)
    {
        // No bonds, add 4 hydrogens in tetrahedral arrangement
        // Use standard tetrahedral vertices
        Hydrogens.Add(Center + FVector(1, 1, 1).GetSafeNormal() * BondLength);
        Hydrogens.Add(Center + FVector(1, -1, -1).GetSafeNormal() * BondLength);
        Hydrogens.Add(Center + FVector(-1, 1, -1).GetSafeNormal() * BondLength);
        Hydrogens.Add(Center + FVector(-1, -1, 1).GetSafeNormal() * BondLength);
    }
    
    return Hydrogens;
}

TArray<FVector> FHydrogenGenerator::GenerateTrigonalHydrogens(
    const FVector& Center,
    const TArray<FVector>& Existing,
    int32 NumH,
    float BondLength)
{
    TArray<FVector> Hydrogens;
    
    // Trigonal planar: 120 degrees
    const float TrigAngle = FMath::DegreesToRadians(120.0f);
    
    if (Existing.Num() == 2 && NumH == 1)
    {
        // 2 bonds exist, add 1 hydrogen in plane
        FVector V1 = (Existing[0] - Center).GetSafeNormal();
        FVector V2 = (Existing[1] - Center).GetSafeNormal();
        
        // Get normal to plane
        FVector Normal = FVector::CrossProduct(V1, V2).GetSafeNormal();
        if (Normal.SizeSquared() < KINDA_SMALL_NUMBER)
            Normal = GetPerpendicularVector(V1);
        
        // Bisector of existing bonds (pointing opposite direction)
        FVector Bisector = -(V1 + V2).GetSafeNormal();
        Hydrogens.Add(Center + Bisector * BondLength);
    }
    else if (Existing.Num() == 1 && NumH == 2)
    {
        // 1 bond exists, add 2 hydrogens at 120 degrees
        FVector V = (Existing[0] - Center).GetSafeNormal();
        FVector Perp = GetPerpendicularVector(V);
        FVector Normal = FVector::CrossProduct(V, Perp).GetSafeNormal();
        
        // Two hydrogens at 120 degrees
        FQuat Q1 = FQuat(Normal, TrigAngle);
        FQuat Q2 = FQuat(Normal, -TrigAngle);
        
        Hydrogens.Add(Center + Q1.RotateVector(V) * BondLength);
        Hydrogens.Add(Center + Q2.RotateVector(V) * BondLength);
    }
    else if (Existing.Num() == 0 && NumH == 3)
    {
        // No bonds, add 3 hydrogens in trigonal planar
        FVector Up = FVector::UpVector;
        FVector Perp = GetPerpendicularVector(Up);
        
        for (int32 i = 0; i < 3; ++i)
        {
            float Angle = (i * 120.0f) * PI / 180.0f;
            FQuat Rotation = FQuat(Up, Angle);
            FVector Dir = Rotation.RotateVector(Perp);
            Hydrogens.Add(Center + Dir * BondLength);
        }
    }
    
    return Hydrogens;
}

TArray<FVector> FHydrogenGenerator::GenerateLinearHydrogens(
    const FVector& Center,
    const TArray<FVector>& Existing,
    int32 NumH,
    float BondLength)
{
    TArray<FVector> Hydrogens;
    
    if (Existing.Num() == 1 && NumH == 1)
    {
        // 1 bond exists, add 1 hydrogen on opposite side (linear)
        FVector V = (Existing[0] - Center).GetSafeNormal();
        Hydrogens.Add(Center - V * BondLength);
    }
    else if (Existing.Num() == 0 && NumH == 2)
    {
        // No bonds, add 2 hydrogens in opposite directions
        FVector Dir = FVector::ForwardVector;
        Hydrogens.Add(Center + Dir * BondLength);
        Hydrogens.Add(Center - Dir * BondLength);
    }
    
    return Hydrogens;
}

FVector FHydrogenGenerator::GetPerpendicularVector(const FVector& V)
{
    FVector Perp = FVector::CrossProduct(V, FVector::UpVector);
    if (Perp.SizeSquared() < KINDA_SMALL_NUMBER)
        Perp = FVector::CrossProduct(V, FVector::RightVector);
    return Perp.GetSafeNormal();
}