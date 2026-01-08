// PDBTreeEntryWidget.cpp
#include "PDBTreeEntryWidget.h"
#include "PDBViewer.h"
#include "Components/TextBlock.h"
#include "Components/Button.h"
#include "Components/CheckBox.h"
#include "Components/TreeView.h"
#include "Components/PanelWidget.h"
#include "Kismet/GameplayStatics.h"

void UPDBTreeEntryWidget::NativeOnListItemObjectSet(UObject* ListItemObject)
{
    IUserObjectListEntry::NativeOnListItemObjectSet(ListItemObject);
    
    // Cast to our node type
    CurrentNode = Cast<UPDBTreeNode>(ListItemObject);
    if (!CurrentNode)
        return;
    
    // Find PDBViewer if we don't have it
    if (!PDBViewerRef)
    {
        TArray<AActor*> FoundActors;
        UGameplayStatics::GetAllActorsOfClass(GetWorld(), APDBViewer::StaticClass(), FoundActors);
        if (FoundActors.Num() > 0)
        {
            PDBViewerRef = Cast<APDBViewer>(FoundActors[0]);
        }
    }
    
    // Set display name with proper color
    if (DisplayNameText)
    {
        DisplayNameText->SetText(FText::FromString(CurrentNode->DisplayName));
        // Auto-set text color to primary
        DisplayNameText->SetColorAndOpacity(FSlateColor(FLinearColor::FromSRGBColor(FColor::FromHex(TEXT("1F2937")))));
    }
    
    // Set visibility checkbox
    if (VisibilityCheckbox)
    {
        VisibilityCheckbox->SetIsChecked(CurrentNode->bIsVisible);
        
        // Style the checkbox for better visibility
        FCheckBoxStyle CheckboxStyle = VisibilityCheckbox->GetWidgetStyle();
        
        // Checked state - bright pink
        CheckboxStyle.CheckedImage.TintColor = FSlateColor(FLinearColor::FromSRGBColor(FColor::FromHex(TEXT("EC4899"))));
        CheckboxStyle.CheckedHoveredImage.TintColor = FSlateColor(FLinearColor::FromSRGBColor(FColor::FromHex(TEXT("DB2777"))));
        CheckboxStyle.CheckedPressedImage.TintColor = FSlateColor(FLinearColor::FromSRGBColor(FColor::FromHex(TEXT("BE185D"))));
        
        // Unchecked state - gray border
        CheckboxStyle.UncheckedImage.TintColor = FSlateColor(FLinearColor::FromSRGBColor(FColor::FromHex(TEXT("9CA3AF"))));
        CheckboxStyle.UncheckedHoveredImage.TintColor = FSlateColor(FLinearColor::FromSRGBColor(FColor::FromHex(TEXT("6B7280"))));
        CheckboxStyle.UncheckedPressedImage.TintColor = FSlateColor(FLinearColor::FromSRGBColor(FColor::FromHex(TEXT("4B5563"))));
        
        VisibilityCheckbox->SetWidgetStyle(CheckboxStyle);
        
        // Bind checkbox event
        if (!VisibilityCheckbox->OnCheckStateChanged.IsBound())
        {
            VisibilityCheckbox->OnCheckStateChanged.AddDynamic(this, &UPDBTreeEntryWidget::OnVisibilityChanged);
        }
    }
    
    // Show/hide expander button based on whether it's a chain
    if (ExpanderButton)
    {
        ExpanderButton->SetVisibility(CurrentNode->bIsChain ? ESlateVisibility::Visible : ESlateVisibility::Collapsed);
        
        // Set expander button text color
        UPanelWidget* ButtonContent = Cast<UPanelWidget>(ExpanderButton->GetChildAt(0));
        if (ButtonContent && ButtonContent->GetChildrenCount() > 0)
        {
            if (UTextBlock* ButtonText = Cast<UTextBlock>(ButtonContent->GetChildAt(0)))
            {
                ButtonText->SetColorAndOpacity(FSlateColor(FLinearColor::FromSRGBColor(FColor::FromHex(TEXT("6B7280")))));
            }
        }
        
        // Bind expander button event
        if (!ExpanderButton->OnClicked.IsBound())
        {
            ExpanderButton->OnClicked.AddDynamic(this, &UPDBTreeEntryWidget::OnExpanderClicked);
        }
    }
}

void UPDBTreeEntryWidget::OnExpanderClicked()
{
    if (!CurrentNode)
        return;
    
    // Get the owning tree view
    UTreeView* TreeView = Cast<UTreeView>(GetOwningListView());
    if (TreeView)
    {
        // Toggle the expansion state
        bIsExpanded = !bIsExpanded;
        TreeView->SetItemExpansion(CurrentNode, bIsExpanded);
        
        // Update button text/icon
        if (ExpanderButton)
        {
            UTextBlock* ButtonText = Cast<UTextBlock>(ExpanderButton->GetChildAt(0));
            if (ButtonText)
            {
                ButtonText->SetText(FText::FromString(bIsExpanded ? TEXT("▼") : TEXT("▶")));
            }
        }
    }
}

void UPDBTreeEntryWidget::OnVisibilityChanged(bool bIsChecked)
{
    if (!CurrentNode || !PDBViewerRef)
        return;
    
    // Toggle visibility in the PDBViewer
    PDBViewerRef->ToggleNodeVisibility(CurrentNode);
}