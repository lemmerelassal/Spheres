// HydrogenGenerator.cpp - Implementation of hydrogen generation

#include "HydrogenGenerator.h"

TArray<FVector> FHydrogenGenerator::GenerateHydrogenPositions(
    const FVector& ParentPos,
    int32 ParentIdx,
    int32 NumHydrogens,
    const TArray<FVector>& AllAtomPositions,
    const TArray<TPair<int32, int32>>& BondPairs,
    const FString& ParentElement)
{
    TArray<FVector> HPositions;
    
    // Use the ideal bond length function instead of hardcoded value
    const float BondLength = GetIdealBondLength(ParentElement);
    
    // Find existing bond directions from this atom
    TArray<FVector> ExistingBondDirs;
    for (const auto& Bond : BondPairs)
    {
        int32 OtherIdx = -1;
        if (Bond.Key == ParentIdx)
            OtherIdx = Bond.Value;
        else if (Bond.Value == ParentIdx)
            OtherIdx = Bond.Key;
        
        if (OtherIdx >= 0 && AllAtomPositions.IsValidIndex(OtherIdx))
        {
            FVector Dir = (AllAtomPositions[OtherIdx] - ParentPos).GetSafeNormal();
            ExistingBondDirs.Add(Dir);
        }
    }
    
    FString Element = ParentElement.ToUpper();
    
    // Generate hydrogen positions based on geometry
    if (NumHydrogens == 1)
    {
        // Single hydrogen - geometry depends on element and existing bonds
        if (ExistingBondDirs.Num() == 1)
        {
            FVector V = ExistingBondDirs[0];
            
            // For oxygen/sulfur: use bent geometry (like water), not linear
            if (Element == TEXT("O") || Element == TEXT("S"))
            {
                // Bent geometry: ~104.5° for oxygen, ~92° for sulfur
                float BentAngle = (Element == TEXT("O")) ? 
                    FMath::DegreesToRadians(104.5f) : FMath::DegreesToRadians(92.0f);
                
                // Get perpendicular to create the bent angle
                FVector Perp = GetPerpendicularVector(V);
                
                // Rotate from -V by half the bent angle
                float HalfAngle = BentAngle / 2.0f;
                FVector HDir = (-V * FMath::Cos(HalfAngle) + Perp * FMath::Sin(HalfAngle)).GetSafeNormal();
                
                HPositions.Add(ParentPos + HDir * BondLength);
            }
            else
            {
                // For other elements (C, N, etc): use tetrahedral-like geometry
                // Place at tetrahedral angle from the existing bond
                FVector Perp = GetPerpendicularVector(V);
                float TetAngle = FMath::DegreesToRadians(109.471f);
                FVector HDir = (V * FMath::Cos(TetAngle) + Perp * FMath::Sin(TetAngle)).GetSafeNormal();
                HPositions.Add(ParentPos + HDir * BondLength);
            }
        }
        else
        {
            // Multiple existing bonds - place opposite to sum
            FVector SumDir = FVector::ZeroVector;
            for (const FVector& Dir : ExistingBondDirs)
                SumDir += Dir;
            
            FVector HDir = -SumDir.GetSafeNormal();
            if (HDir.SizeSquared() < 0.01f) // If no existing bonds or they cancel
                HDir = FVector::UpVector;
            
            HPositions.Add(ParentPos + HDir * BondLength);
        }
    }
    else if (NumHydrogens == 2)
    {
        // Two hydrogens - geometry depends on element and existing bonds
        if (ExistingBondDirs.Num() == 2)
        {
            // Two existing bonds + 2 hydrogens (rare case)
            if (Element == TEXT("O") || Element == TEXT("S"))
            {
                // Water-like with 2 other bonds - use tetrahedral bent
                FVector SumDir = (ExistingBondDirs[0] + ExistingBondDirs[1]).GetSafeNormal();
                FVector Perp = FVector::CrossProduct(ExistingBondDirs[0], ExistingBondDirs[1]).GetSafeNormal();
                
                if (Perp.SizeSquared() < 0.01f)
                    Perp = FVector::UpVector;
                
                // Use ~104.5° for oxygen (half angle = 52.25°)
                float HalfAngle = (Element == TEXT("O")) ? 52.25f : 60.0f;
                FVector Base = -SumDir.GetSafeNormal();
                
                HPositions.Add(ParentPos + (Base.RotateAngleAxis(HalfAngle, Perp)).GetSafeNormal() * BondLength);
                HPositions.Add(ParentPos + (Base.RotateAngleAxis(-HalfAngle, Perp)).GetSafeNormal() * BondLength);
            }
            else if (Element == TEXT("N"))
            {
                // NH2 with 2 existing bonds - use bent geometry (~107°)
                FVector SumDir = (ExistingBondDirs[0] + ExistingBondDirs[1]).GetSafeNormal();
                FVector Perp = FVector::CrossProduct(ExistingBondDirs[0], ExistingBondDirs[1]).GetSafeNormal();
                
                if (Perp.SizeSquared() < 0.01f)
                    Perp = FVector::UpVector;
                
                FVector Base = -SumDir.GetSafeNormal();
                
                HPositions.Add(ParentPos + (Base.RotateAngleAxis(53.5f, Perp)).GetSafeNormal() * BondLength);
                HPositions.Add(ParentPos + (Base.RotateAngleAxis(-53.5f, Perp)).GetSafeNormal() * BondLength);
            }
            else
            {
                // Carbon or other - tetrahedral
                FVector SumDir = (ExistingBondDirs[0] + ExistingBondDirs[1]).GetSafeNormal();
                FVector Perp = FVector::CrossProduct(ExistingBondDirs[0], ExistingBondDirs[1]).GetSafeNormal();
                
                if (Perp.SizeSquared() < 0.01f)
                    Perp = FVector::UpVector;
                
                FVector Base = -SumDir.GetSafeNormal();
                
                HPositions.Add(ParentPos + (Base.RotateAngleAxis(60.0f, Perp)).GetSafeNormal() * BondLength);
                HPositions.Add(ParentPos + (Base.RotateAngleAxis(-60.0f, Perp)).GetSafeNormal() * BondLength);
            }
        }
        else if (ExistingBondDirs.Num() == 1)
        {
            // One existing bond + 2 hydrogens (most common case)
            FVector V = ExistingBondDirs[0];
            FVector Perp = GetPerpendicularVector(V);
            
            if (Element == TEXT("O") || Element == TEXT("S"))
            {
                // Water-like (OH2, SH2) - use bent geometry with proper angle
                // For water: H-O-H angle is ~104.5°
                // We need to place hydrogens at 104.5° from each other
                
                FVector Normal = FVector::CrossProduct(V, Perp).GetSafeNormal();
                
                // The angle from the existing bond direction to each H
                // If H-O-H is 104.5°, and we're placing them symmetrically,
                // each H makes an angle with the opposite of the existing bond
                float HalfHOHAngle = FMath::DegreesToRadians((Element == TEXT("O")) ? 104.5f / 2.0f : 92.0f / 2.0f);
                
                // Base direction opposite to existing bond
                FVector Base = -V;
                
                // Create rotation quaternions
                FQuat Q1 = FQuat(Normal, HalfHOHAngle);
                FQuat Q2 = FQuat(Normal, -HalfHOHAngle);
                
                HPositions.Add(ParentPos + Q1.RotateVector(Base) * BondLength);
                HPositions.Add(ParentPos + Q2.RotateVector(Base) * BondLength);
            }
            else if (Element == TEXT("N"))
            {
                // NH2 group - use trigonal planar geometry (120° angles)
                FVector Normal = FVector::CrossProduct(V, Perp).GetSafeNormal();
                
                // Two hydrogens at 120 degrees from the existing bond
                float TrigAngle = FMath::DegreesToRadians(120.0f);
                FQuat Q1 = FQuat(Normal, TrigAngle);
                FQuat Q2 = FQuat(Normal, -TrigAngle);
                
                HPositions.Add(ParentPos + Q1.RotateVector(V) * BondLength);
                HPositions.Add(ParentPos + Q2.RotateVector(V) * BondLength);
            }
            else
            {
                // Carbon or other - use tetrahedral arrangement
                // This would be something like CH2 in a ring or chain
                // For 1 bond + 2 H: need to place 2 H's at tetrahedral angles from the bond
                
                float TetAngle = FMath::DegreesToRadians(109.471f);
                
                // Get perpendicular to the existing bond
                FVector Perp1 = GetPerpendicularVector(V);
                FVector Perp2 = FVector::CrossProduct(V, Perp1).GetSafeNormal();
                
                // For tetrahedral geometry, when you have 1 bond and need to add 2 hydrogens,
                // they should be positioned such that all angles are 109.5°
                // The two hydrogens form a "V" shape that's tilted relative to the existing bond
                
                // Place them using a combination of the two perpendicular directions
                // This creates the proper 3D tetrahedral arrangement
                float CosTheta = FMath::Cos(TetAngle);
                float SinTheta = FMath::Sin(TetAngle);
                
                // The two hydrogens are at ±109.5° from V, symmetric about a plane
                // H1 and H2 are rotated 120° apart around V axis (for proper tetrahedral spacing)
                FVector HDir1 = (V * CosTheta + Perp1 * SinTheta).GetSafeNormal();
                
                // Rotate Perp1 by 120° around V to get the second hydrogen direction
                FQuat RotQuat = FQuat(V, FMath::DegreesToRadians(120.0f));
                FVector Perp1Rotated = RotQuat.RotateVector(Perp1);
                FVector HDir2 = (V * CosTheta + Perp1Rotated * SinTheta).GetSafeNormal();
                
                HPositions.Add(ParentPos + HDir1 * BondLength);
                HPositions.Add(ParentPos + HDir2 * BondLength);
            }
        }
        else
        {
            // No existing bonds - default: opposite directions
            HPositions.Add(ParentPos + FVector::UpVector * BondLength);
            HPositions.Add(ParentPos + FVector::DownVector * BondLength);
        }
    }
    else if (NumHydrogens == 3)
    {
        // Three hydrogens
        if (Element == TEXT("N") && ExistingBondDirs.Num() == 1)
        {
            // NH3+ (ammonium) or similar - use tetrahedral
            FVector V = ExistingBondDirs[0];
            FVector Perp = GetPerpendicularVector(V);
            
            // Rotate 120 degrees around the existing bond
            for (int32 i = 0; i < 3; ++i)
            {
                float Angle = (i * 120.0f) * PI / 180.0f;
                FQuat Rotation = FQuat(V, Angle);
                FVector Rotated = Rotation.RotateVector(Perp);
                
                // Tilt down by tetrahedral angle
                float TetAngle = FMath::DegreesToRadians(109.471f);
                FVector Cross = FVector::CrossProduct(V, Rotated).GetSafeNormal();
                FQuat Tilt = FQuat(Cross, PI - TetAngle);
                FVector Dir = Tilt.RotateVector(V);
                
                HPositions.Add(ParentPos + Dir * BondLength);
            }
        }
        else
        {
            // CH3 or similar - tetrahedral geometry
            if (ExistingBondDirs.Num() > 0)
            {
                // V points from parent atom to bonded atom (away from CH3)
                FVector V = ExistingBondDirs[0];
                FVector Perp = GetPerpendicularVector(V);
                
                // The hydrogens should point in same general direction as V, tilted by tetrahedral angle
                // This creates a tetrahedral arrangement where H's are on same side as the bond
                float TetAngle = FMath::DegreesToRadians(109.471f);
                
                // Rotate perpendicular vector around V axis to get 3 positions at 120° intervals
                for (int32 i = 0; i < 3; ++i)
                {
                    float Angle = (i * 120.0f) * PI / 180.0f;
                    FQuat Rotation = FQuat(V, Angle);
                    FVector RotatedPerp = Rotation.RotateVector(Perp);
                    
                    // Tilt from V by tetrahedral angle toward the rotated perpendicular direction
                    FVector HDir = (V * FMath::Cos(TetAngle) + RotatedPerp * FMath::Sin(TetAngle)).GetSafeNormal();
                    
                    HPositions.Add(ParentPos + HDir * BondLength);
                }
            }
            else
            {
                // No existing bond - use default tetrahedral arrangement
                FVector BaseDir = FVector::UpVector;
                FVector Perp = GetPerpendicularVector(BaseDir);
                
                float TetAngle = FMath::DegreesToRadians(109.471f);
                
                for (int32 i = 0; i < 3; ++i)
                {
                    float Angle = (i * 120.0f) * PI / 180.0f;
                    FQuat Rotation = FQuat(BaseDir, Angle);
                    FVector RotatedPerp = Rotation.RotateVector(Perp);
                    
                    FVector HDir = (BaseDir * FMath::Cos(TetAngle) + RotatedPerp * FMath::Sin(TetAngle)).GetSafeNormal();
                    
                    HPositions.Add(ParentPos + HDir * BondLength);
                }
            }
        }
    }
    else if (NumHydrogens == 4)
    {
        // Four hydrogens - methane tetrahedral geometry
        HPositions.Add(ParentPos + FVector(1, 1, 1).GetSafeNormal() * BondLength);
        HPositions.Add(ParentPos + FVector(-1, -1, 1).GetSafeNormal() * BondLength);
        HPositions.Add(ParentPos + FVector(-1, 1, -1).GetSafeNormal() * BondLength);
        HPositions.Add(ParentPos + FVector(1, -1, -1).GetSafeNormal() * BondLength);
    }
    
    return HPositions;
}

// Helper function to get a perpendicular vector
FVector FHydrogenGenerator::GetPerpendicularVector(const FVector& V)
{
    FVector Perp = FVector::CrossProduct(V, FVector::UpVector);
    if (Perp.SizeSquared() < KINDA_SMALL_NUMBER)
        Perp = FVector::CrossProduct(V, FVector::RightVector);
    return Perp.GetSafeNormal();
}

int32 FHydrogenGenerator::GetTypicalValence(const FString &Element)
{
    static const TMap<FString, int32> Valences = {
        {TEXT("C"), 4},  // Carbon
        {TEXT("N"), 3},  // Nitrogen (can be 3 or 5, use 3 for organic)
        {TEXT("O"), 2},  // Oxygen
        {TEXT("S"), 2},  // Sulfur (can be 2, 4, or 6, use 2 for organic)
        {TEXT("P"), 3},  // Phosphorus (can be 3 or 5)
        {TEXT("F"), 1},  // Fluorine
        {TEXT("CL"), 1}, // Chlorine
        {TEXT("BR"), 1}, // Bromine
        {TEXT("I"), 1},  // Iodine
        {TEXT("B"), 3},  // Boron
        {TEXT("SI"), 4}, // Silicon
    };

    const int32 *ValencePtr = Valences.Find(Element);
    return ValencePtr ? *ValencePtr : 0;
}

TArray<TPair<FVector, int32>> FHydrogenGenerator::GenerateHydrogens(
    const TArray<FVector> &AtomPositions,
    const TArray<FString> &AtomElements,
    const TArray<TPair<int32, int32>> &BondPairs,
    const TArray<int32> &BondOrders)
{
    TArray<TPair<FVector, int32>> Hydrogens;

    if (AtomPositions.Num() != AtomElements.Num())
    {
        UE_LOG(LogTemp, Error, TEXT("HydrogenGenerator: Atom position/element count mismatch"));
        return Hydrogens;
    }

    // Count bonds for each atom
    TArray<int32> BondCounts;
    BondCounts.SetNum(AtomPositions.Num());
    for (int32 i = 0; i < BondCounts.Num(); ++i)
        BondCounts[i] = 0;

    // Count total bond order for each atom
    for (int32 i = 0; i < BondPairs.Num(); ++i)
    {
        int32 Atom1 = BondPairs[i].Key;
        int32 Atom2 = BondPairs[i].Value;
        int32 Order = BondOrders.IsValidIndex(i) ? BondOrders[i] : 1;

        // Treat aromatic bonds (order 4) as 1.5 bonds, round up to 2
        if (Order == 4)
            Order = 2;

        if (BondCounts.IsValidIndex(Atom1))
            BondCounts[Atom1] += Order;
        if (BondCounts.IsValidIndex(Atom2))
            BondCounts[Atom2] += Order;
    }

    UE_LOG(LogTemp, Log, TEXT("HydrogenGenerator: Processing %d atoms, %d bonds"),
           AtomPositions.Num(), BondPairs.Num());

    // Generate hydrogens for each heavy atom
    for (int32 AtomIdx = 0; AtomIdx < AtomPositions.Num(); ++AtomIdx)
    {
        FString Element = AtomElements[AtomIdx].ToUpper();

        // Skip if already hydrogen
        if (Element == TEXT("H") || Element == TEXT("D"))
            continue;

        // Get typical valence for this element
        int32 TypicalValence = GetTypicalValence(Element);
        if (TypicalValence == 0)
            continue; // Unknown element or doesn't typically bond to H

        // Calculate how many hydrogens needed
        int32 CurrentBonds = BondCounts[AtomIdx];
        int32 HydrogensNeeded = TypicalValence - CurrentBonds;

        // Log detailed info for each atom
        UE_LOG(LogTemp, Log, TEXT("  Atom %d (%s): Valence=%d, CurrentBonds=%d, HNeeded=%d"),
               AtomIdx, *Element, TypicalValence, CurrentBonds, HydrogensNeeded);

        if (HydrogensNeeded <= 0)
            continue;

        // Clamp to reasonable values
        HydrogensNeeded = FMath::Clamp(HydrogensNeeded, 0, 4);

        // Generate hydrogen positions
        TArray<FVector> HPositions = GenerateHydrogenPositions(
            AtomPositions[AtomIdx],
            AtomIdx,
            HydrogensNeeded,
            AtomPositions,
            BondPairs,
            Element);

        // Add each hydrogen
        for (const FVector &HPos : HPositions)
        {
            Hydrogens.Add(TPair<FVector, int32>(HPos, AtomIdx));
        }

        UE_LOG(LogTemp, Log, TEXT("    Generated %d hydrogens for atom %d"),
               HPositions.Num(), AtomIdx);
    }

    UE_LOG(LogTemp, Log, TEXT("HydrogenGenerator: Total hydrogens generated: %d"), Hydrogens.Num());

    return Hydrogens;
}

int32 FHydrogenGenerator::GetIdealHydrogenCount(const FString &Element, int32 TotalBondOrder)
{
    FString Elem = Element.ToUpper();

    // Calculate based on total bond order (sum of all bond orders to this atom)
    // Valence - TotalBondOrder = number of H needed

    if (Elem == TEXT("C"))
        return FMath::Max(0, 4 - TotalBondOrder); // sp³: 4, sp²: 3, sp: 2
    else if (Elem == TEXT("N"))
        return FMath::Max(0, 3 - TotalBondOrder); // Can have 0-3 H depending on bonds
    else if (Elem == TEXT("O"))
        return FMath::Max(0, 2 - TotalBondOrder); // Can have 0-2 H
    else if (Elem == TEXT("S"))
        return FMath::Max(0, 2 - TotalBondOrder); // Similar to oxygen
    else if (Elem == TEXT("P"))
        return FMath::Max(0, 3 - TotalBondOrder); // Can have 0-3 H
    else if (Elem == TEXT("F") || Elem == TEXT("CL") ||
             Elem == TEXT("BR") || Elem == TEXT("I"))
        return FMath::Max(0, 1 - TotalBondOrder); // Halogens

    return 0;
}

float FHydrogenGenerator::GetIdealBondLength(const FString &Element)
{
    FString Elem = Element.ToUpper();

    // Bond lengths in Angstroms, scaled by PDB::SCALE (50.0)
    if (Elem == TEXT("C"))
        return 1.09f; // C-H
    else if (Elem == TEXT("N"))
        return 1.01f; // N-H
    else if (Elem == TEXT("O"))
        return 0.96f; // O-H
    else if (Elem == TEXT("S"))
        return 1.34f; // S-H
    else if (Elem == TEXT("P"))
        return 1.42f; // P-H

    return 1.09f; // Default
}

TArray<FVector> FHydrogenGenerator::GenerateHydrogensForAtom(
    const FVector &AtomPos,
    const FString &Element,
    const TArray<FVector> &BondedPositions,
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
    const FVector &Center,
    const TArray<FVector> &Existing,
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

        // Bisector of existing bonds (average direction)
        FVector Bisector = (V1 + V2).GetSafeNormal();
        
        // For tetrahedral geometry, we want the hydrogens to be placed such that
        // all 4 bond angles are 109.5°. The hydrogens should be symmetric
        // about the plane containing the two existing bonds.
        
        // Create a vector perpendicular to the bisector, lying in the perpendicular plane
        FVector PerpendicularInPlane = FVector::CrossProduct(Normal, Bisector).GetSafeNormal();
        
        // The angle between the existing bond bisector and the H direction
        // For perfect tetrahedral: if two bonds are at angle θ apart,
        // the hydrogens should complete the tetrahedron
        float ExistingAngle = FMath::Acos(FMath::Clamp(FVector::DotProduct(V1, V2), -1.0f, 1.0f));
        
        // For tetrahedral, the angle from -Bisector to each H should be calculated
        // to maintain 109.5° with the existing bonds
        // Using the formula: for tetrahedral geometry with 2+2 arrangement
        float HAngleFromBisector = FMath::Acos((-FMath::Cos(ExistingAngle) - FMath::Cos(TetAngle)) / 
                                                 (1 + FMath::Cos(ExistingAngle)));
        
        // Direction opposite to bisector, then spread the hydrogens
        FVector BaseDir = -Bisector;
        
        // Rotate BaseDir toward the perpendicular by the calculated angle
        FQuat Q1 = FQuat(PerpendicularInPlane, HAngleFromBisector);
        FQuat Q2 = FQuat(PerpendicularInPlane, -HAngleFromBisector);

        Hydrogens.Add(Center + Q1.RotateVector(BaseDir) * BondLength);
        Hydrogens.Add(Center + Q2.RotateVector(BaseDir) * BondLength);
    }
    else if (Existing.Num() == 1 && NumH == 3)
    {
        // 1 bond exists, add 3 hydrogens
        // V points from center to the existing bonded atom
        FVector V = (Existing[0] - Center).GetSafeNormal();
        FVector Perp = GetPerpendicularVector(V);

        // The hydrogens should point in same general direction as V, tilted by tetrahedral angle
        // Rotate perpendicular vector around V axis to get 3 positions at 120° intervals
        for (int32 i = 0; i < 3; ++i)
        {
            float Angle = (i * 120.0f) * PI / 180.0f;
            FQuat Rotation = FQuat(V, Angle);
            FVector RotatedPerp = Rotation.RotateVector(Perp);
            
            // Tilt from V by tetrahedral angle
            FVector HDir = (V * FMath::Cos(TetAngle) + RotatedPerp * FMath::Sin(TetAngle)).GetSafeNormal();
            
            Hydrogens.Add(Center + HDir * BondLength);
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
    const FVector &Center,
    const TArray<FVector> &Existing,
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
    const FVector &Center,
    const TArray<FVector> &Existing,
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
