// MoleculeListEntry.cpp
#include "MoleculeListEntry.h"
#include "PDBViewer.h"
#include "MMGBSA.h"
#include "Components/TextBlock.h"
#include "Components/Button.h"
#include "Kismet/GameplayStatics.h"

void UMoleculeListEntry::NativeConstruct()
{
    Super::NativeConstruct();

    if (ButtonToggleVisibility)
    {
        ButtonToggleVisibility->OnClicked.AddDynamic(this, &UMoleculeListEntry::OnToggleClicked);
    }
    
    if (ButtonCalculateAffinity)
    {
        ButtonCalculateAffinity->OnClicked.AddDynamic(this, &UMoleculeListEntry::OnCalculateAffinityClicked);
    }
    
    // Find or create MMGBSA actor
    TArray<AActor*> Found;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), AMMGBSA::StaticClass(), Found);
    if (Found.Num() > 0)
    {
        MMGBSARef = Cast<AMMGBSA>(Found[0]);
    }
    else
    {
        // Spawn MMGBSA actor
        MMGBSARef = GetWorld()->SpawnActor<AMMGBSA>(AMMGBSA::StaticClass());
        
        // Initialize from viewer
        APDBViewer* Viewer = Cast<APDBViewer>(
            UGameplayStatics::GetActorOfClass(GetWorld(), APDBViewer::StaticClass())
        );
        if (Viewer && MMGBSARef)
        {
            MMGBSARef->InitializeFromViewer(Viewer);
        }
    }
    
    // Bind to affinity calculation event
    if (MMGBSARef)
    {
        MMGBSARef->OnAffinityCalculated.AddDynamic(this, &UMoleculeListEntry::OnAffinityCalculated);
    }

    ApplyStyling();
}

void UMoleculeListEntry::ApplyStyling()
{
    // Style the molecule name
    if (TextMoleculeName)
    {
        TextMoleculeName->SetColorAndOpacity(FLinearColor(0.9f, 0.9f, 0.9f));
    }

    // Style the atom and bond counts
    if (TextAtomCount)
    {
        TextAtomCount->SetColorAndOpacity(FLinearColor(0.7f, 0.7f, 0.7f));
    }

    if (TextBondCount)
    {
        TextBondCount->SetColorAndOpacity(FLinearColor(0.7f, 0.7f, 0.7f));
    }

    // Style the toggle button
    if (ButtonToggleVisibility)
    {
        ButtonToggleVisibility->SetColorAndOpacity(FLinearColor(0.15f, 0.3f, 0.6f));
        
        // Style the button text
        UTextBlock* ButtonText = Cast<UTextBlock>(ButtonToggleVisibility->GetChildAt(0));
        if (ButtonText)
        {
            ButtonText->SetColorAndOpacity(FLinearColor::White);
        }
    }
}

void UMoleculeListEntry::NativeOnListItemObjectSet(UObject* ListItemObject)
{
    // This is called when the ListView assigns a molecule node to this entry
    CurrentMoleculeNode = Cast<UPDBMoleculeNode>(ListItemObject);

    if (!CurrentMoleculeNode)
        return;

    // Update the UI elements
    if (TextMoleculeName)
    {
        TextMoleculeName->SetText(FText::FromString(CurrentMoleculeNode->MoleculeName));
    }

    if (TextAtomCount)
    {
        FString AtomText = FString::Printf(TEXT("Atoms: %d"), CurrentMoleculeNode->AtomCount);
        TextAtomCount->SetText(FText::FromString(AtomText));
    }

    if (TextBondCount)
    {
        FString BondText = FString::Printf(TEXT("Bonds: %d"), CurrentMoleculeNode->BondCount);
        TextBondCount->SetText(FText::FromString(BondText));
    }

    // Update button appearance
    UpdateButtonText();
    
    // Check for cached affinity result
    UpdateAffinityDisplay();
}

void UMoleculeListEntry::OnToggleClicked()
{
    if (!CurrentMoleculeNode)
        return;

    APDBViewer* Viewer = Cast<APDBViewer>(
        UGameplayStatics::GetActorOfClass(GetWorld(), APDBViewer::StaticClass())
    );

    if (Viewer)
    {
        Viewer->ToggleMoleculeNodeVisibility(CurrentMoleculeNode);
        UpdateButtonText();
    }
}

void UMoleculeListEntry::OnCalculateAffinityClicked()
{
    if (!CurrentMoleculeNode || !MMGBSARef)
        return;
    
    if (TextBindingAffinity)
    {
        TextBindingAffinity->SetText(FText::FromString(TEXT("Calculating...")));
        TextBindingAffinity->SetColorAndOpacity(FSlateColor(FLinearColor::Gray));
    }

    // Calculate binding affinity for this ligand
    FBindingAffinityResult Result = MMGBSARef->CalculateBindingAffinity(CurrentMoleculeNode->MoleculeKey);
    
    if (Result.bIsValid)
    {
        UpdateAffinityDisplay();
    }
    else
    {
        if (TextBindingAffinity)
        {
            TextBindingAffinity->SetText(FText::FromString(TEXT("Calculation failed")));
            TextBindingAffinity->SetColorAndOpacity(FSlateColor(FLinearColor::Red));
        }
        UE_LOG(LogTemp, Warning, TEXT("MoleculeListEntry: MMGBSA failed for %s"), *CurrentMoleculeNode->MoleculeKey);
    }
}

void UMoleculeListEntry::OnAffinityCalculated(const FBindingAffinityResult& Result)
{
    // Check if this result is for our molecule
    if (CurrentMoleculeNode && Result.LigandKey == CurrentMoleculeNode->MoleculeKey)
    {
        UpdateAffinityDisplay();
    }
}

void UMoleculeListEntry::UpdateButtonText()
{
    if (!CurrentMoleculeNode || !ButtonToggleVisibility)
        return;

    UTextBlock* ButtonText = Cast<UTextBlock>(ButtonToggleVisibility->GetChildAt(0));
    if (ButtonText)
    {
        FString ToggleText = CurrentMoleculeNode->bIsVisible ? TEXT("Hide") : TEXT("Show");
        ButtonText->SetText(FText::FromString(ToggleText));

        // Change button color based on visibility state
        if (CurrentMoleculeNode->bIsVisible)
        {
            ButtonToggleVisibility->SetColorAndOpacity(FLinearColor(0.2f, 0.5f, 0.2f)); // Green when visible
        }
        else
        {
            ButtonToggleVisibility->SetColorAndOpacity(FLinearColor(0.5f, 0.2f, 0.2f)); // Red when hidden
        }
    }
}

void UMoleculeListEntry::UpdateAffinityDisplay()
{
    if (!CurrentMoleculeNode || !MMGBSARef)
        return;
    
    // Check for cached result
    FBindingAffinityResult Result;
    if (MMGBSARef->GetCachedResult(CurrentMoleculeNode->MoleculeKey, Result))
    {
        if (!Result.bIsValid)
        {
            // Calculation was attempted but result marked invalid
            if (TextBindingAffinity)
            {
                TextBindingAffinity->SetText(FText::FromString(TEXT("Calculation failed")));
                TextBindingAffinity->SetColorAndOpacity(FSlateColor(FLinearColor::Red));
            }
            if (TextKi) TextKi->SetText(FText::FromString(TEXT("")));
            if (TextAffinityClass) TextAffinityClass->SetText(FText::FromString(TEXT("")));
            return;
        }

        // Display ΔG
        if (TextBindingAffinity)
        {
            FString AffinityText = FString::Printf(TEXT("ΔG: %.2f kcal/mol"), Result.DeltaG_Binding);
            TextBindingAffinity->SetText(FText::FromString(AffinityText));
            
            // Color code by affinity
            FLinearColor Color = FLinearColor::White;
            if (Result.DeltaG_Binding < -10.0f)
                Color = FLinearColor::Green; // Very strong
            else if (Result.DeltaG_Binding < -7.0f)
                Color = FLinearColor(0.5f, 1.0f, 0.5f); // Strong
            else if (Result.DeltaG_Binding < -5.0f)
                Color = FLinearColor::Yellow; // Moderate
            else if (Result.DeltaG_Binding < -3.0f)
                Color = FLinearColor(1.0f, 0.5f, 0.0f); // Weak
            else
                Color = FLinearColor::Red; // Very weak
            
            TextBindingAffinity->SetColorAndOpacity(FSlateColor(Color));
        }
        
        // Display Ki
        if (TextKi)
        {
            FString KiText;
            if (Result.Ki_uM < 0.001f)
                KiText = FString::Printf(TEXT("Ki: %.2f pM"), Result.Ki_uM * 1000000.0f);
            else if (Result.Ki_uM < 1.0f)
                KiText = FString::Printf(TEXT("Ki: %.2f nM"), Result.Ki_uM * 1000.0f);
            else if (Result.Ki_uM < 1000.0f)
                KiText = FString::Printf(TEXT("Ki: %.2f µM"), Result.Ki_uM);
            else
                KiText = FString::Printf(TEXT("Ki: %.2f mM"), Result.Ki_uM / 1000.0f);
            
            TextKi->SetText(FText::FromString(KiText));
        }
        
        // Display affinity class
        if (TextAffinityClass)
        {
            TextAffinityClass->SetText(FText::FromString(Result.AffinityClass));
        }
    }
    else
    {
        // No cached result - show calculation prompt
        if (TextBindingAffinity)
        {
            TextBindingAffinity->SetText(FText::FromString(TEXT("Click to calculate")));
            TextBindingAffinity->SetColorAndOpacity(FSlateColor(FLinearColor::Gray));
        }
        
        if (TextKi)
        {
            TextKi->SetText(FText::FromString(TEXT("")));
        }
        
        if (TextAffinityClass)
        {
            TextAffinityClass->SetText(FText::FromString(TEXT("")));
        }
    }
}