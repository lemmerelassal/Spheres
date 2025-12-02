#pragma once

#include "Widgets/SCompoundWidget.h"
#include "Widgets/DeclarativeSyntaxSupport.h"

DECLARE_DELEGATE_OneParam(FLigandSelectedDelegate, FString);

class SLigandSidebar : public SCompoundWidget
{
public:
    SLATE_BEGIN_ARGS(SLigandSidebar) {}
        SLATE_ARGUMENT(TArray<FString>, LigandIDs)
        SLATE_EVENT(FLigandSelectedDelegate, OnLigandSelected)
    SLATE_END_ARGS()

    void Construct(const FArguments& InArgs);

private:
    TArray<FString> LigandIDs;
    FLigandSelectedDelegate OnLigandSelected;

    TSharedPtr<SListView<TSharedPtr<FString>>> LigandListView;
    TArray<TSharedPtr<FString>> LigandItems;

    TSharedRef<ITableRow> OnGenerateRow(TSharedPtr<FString> Item, const TSharedRef<STableViewBase>& OwnerTable);
    void OnItemClicked(TSharedPtr<FString> Item);
};
