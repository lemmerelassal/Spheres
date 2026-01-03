// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Blueprint/UserWidget.h"
#include "Components/VerticalBox.h"
#include "MyMainWidget.generated.h"

/**
 * 
 */
UCLASS()
class SPHERES_API UMyMainWidget : public UUserWidget
{
	GENERATED_BODY()

	
public:
	// Container that will hold our buttons
	UPROPERTY(meta = (BindWidget))
	UVerticalBox* ButtonContainer;

	// Reference to the button widget Blueprint class
	UPROPERTY(EditDefaultsOnly, BlueprintReadWrite, Category = "UI")
	TSubclassOf<UUserWidget> ButtonWidgetClass;

	// Function to generate buttons
	UFUNCTION(BlueprintCallable, Category = "UI")
	void GenerateButtons(int32 Count);

	// Event that fires when each button is created (implement in Blueprint)
	UFUNCTION(BlueprintImplementableEvent, Category = "UI")
	void OnButtonCreated(UUserWidget* Button, int32 Index);

};
