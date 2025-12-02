#include "SLigandSidebar.h"
#include "Widgets/Input/SButton.h"
#include "Widgets/Text/STextBlock.h"
#include "Widgets/Layout/SScrollBox.h"
#include "Widgets/Views/SListView.h"
#include "Widgets/Layout/SBorder.h"
#include "Widgets/Layout/SBox.h"

void SLigandSidebar::Construct(const FArguments& InArgs)
{
    LigandIDs = InArgs._LigandIDs;
    OnLigandSelected = InArgs._OnLigandSelected;

    for (const FString& ID : LigandIDs)
    {
        LigandItems.Add(MakeShared<FString>(ID));
    }

    ChildSlot
    [
        SNew(SBorder)
        .Padding(8)
        .BorderBackgroundColor(FLinearColor(0.f, 0.f, 0.f, 0.4f))
        [
            SNew(SScrollBox)
            + SScrollBox::Slot()
            [
                SAssignNew(LigandListView, SListView<TSharedPtr<FString>>)
                .ListItemsSource(&LigandItems)
                .OnGenerateRow(this, &SLigandSidebar::OnGenerateRow)
                .SelectionMode(ESelectionMode::Single)
            ]
        ]
    ];
}

TSharedRef<ITableRow> SLigandSidebar::OnGenerateRow(TSharedPtr<FString> Item, const TSharedRef<STableViewBase>& OwnerTable)
{
    return SNew(STableRow<TSharedPtr<FString>>, OwnerTable)
    [
        SNew(SButton)
        .OnClicked_Lambda([this, Item]()
        {
            OnItemClicked(Item);
            return FReply::Handled();
        })
        [
            SNew(STextBlock)
            .Text(FText::FromString(*Item))
            .Justification(ETextJustify::Center)
            .ColorAndOpacity(FSlateColor(FLinearColor::White))
        ]
    ];
}

void SLigandSidebar::OnItemClicked(TSharedPtr<FString> Item)
{
    if (OnLigandSelected.IsBound())
    {
        OnLigandSelected.Execute(*Item);
    }
}
