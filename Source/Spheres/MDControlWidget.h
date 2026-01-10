// MDControlWidget.h - UI for Molecular Dynamics Control
#pragma once

#include "CoreMinimal.h"
#include "Blueprint/UserWidget.h"
#include "MDControlWidget.generated.h"

class AMolecularDynamics;
class UButton;
class UTextBlock;
class USlider;
class UCheckBox;
class UComboBoxString;

UCLASS()
class SPHERES_API UMDControlWidget : public UUserWidget
{
    GENERATED_BODY()

public:
    virtual void NativeConstruct() override;
    virtual void NativeTick(const FGeometry& MyGeometry, float InDeltaTime) override;

protected:
    // Buttons
    UPROPERTY(meta = (BindWidget))
    UButton* Button_Start;
    
    UPROPERTY(meta = (BindWidget))
    UButton* Button_Stop;
    
    UPROPERTY(meta = (BindWidget))
    UButton* Button_Reset;
    
    UPROPERTY(meta = (BindWidget))
    UButton* Button_Initialize;
    
    // Sliders
    UPROPERTY(meta = (BindWidget))
    USlider* Slider_TimeStep;
    
    UPROPERTY(meta = (BindWidget))
    USlider* Slider_Temperature;
    
    UPROPERTY(meta = (BindWidget))
    USlider* Slider_Damping;
    
    // Checkboxes
    UPROPERTY(meta = (BindWidget))
    UCheckBox* CheckBox_LennardJones;
    
    UPROPERTY(meta = (BindWidget))
    UCheckBox* CheckBox_Electrostatics;
    
    UPROPERTY(meta = (BindWidget))
    UCheckBox* CheckBox_Constraints;
    
    // Combo box for integrator selection
    UPROPERTY(meta = (BindWidget))
    UComboBoxString* ComboBox_Integrator;
    
    // Info displays
    UPROPERTY(meta = (BindWidget))
    UTextBlock* Text_Status;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* Text_AtomCount;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* Text_Energy;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* Text_Temperature;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* Text_TimeStep;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* Text_TargetTemp;
    
    UPROPERTY(meta = (BindWidget))
    UTextBlock* Text_DampingValue;
    
    // Binding affinity display (optional)
    UPROPERTY(meta = (BindWidgetOptional))
    UTextBlock* Text_BindingAffinity;

private:
    UPROPERTY()
    AMolecularDynamics* MDSimulation;
    
    // Button callbacks
    UFUNCTION()
    void OnStartClicked();
    
    UFUNCTION()
    void OnStopClicked();
    
    UFUNCTION()
    void OnResetClicked();
    
    UFUNCTION()
    void OnInitializeClicked();
    
    // Slider callbacks
    UFUNCTION()
    void OnTimeStepChanged(float Value);
    
    UFUNCTION()
    void OnTemperatureChanged(float Value);
    
    UFUNCTION()
    void OnDampingChanged(float Value);
    
    // Checkbox callbacks
    UFUNCTION()
    void OnLennardJonesChanged(bool bIsChecked);
    
    UFUNCTION()
    void OnElectrostaticsChanged(bool bIsChecked);
    
    UFUNCTION()
    void OnConstraintsChanged(bool bIsChecked);
    
    // Combobox callback
    UFUNCTION()
    void OnIntegratorChanged(FString SelectedItem, ESelectInfo::Type SelectionType);
    
    // Update displays
    void UpdateInfoDisplays();
    void FindMDSimulation();
};