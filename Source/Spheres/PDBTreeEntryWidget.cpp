// PDBTreeEntryWidget.cpp
#include "PDBTreeEntryWidget.h"
#include "PDBViewer.h"
#include "Components/TextBlock.h"
#include "Components/Button.h"
#include "Components/CheckBox.h"
#include "Components/TreeView.h"
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
    
    // Set display name
    if (DisplayNameText)
    {
        DisplayNameText->SetText(FText::FromString(CurrentNode->DisplayName));
    }
    
    // Set visibility checkbox
    if (VisibilityCheckbox)
    {
        VisibilityCheckbox->SetIsChecked(CurrentNode->bIsVisible);
        
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