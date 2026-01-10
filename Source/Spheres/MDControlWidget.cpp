// MDControlWidget.cpp - UI Implementation for MD Control
#include "MDControlWidget.h"
#include "MolecularDynamics.h"
#include "PDBViewer.h"
#include "Components/Button.h"
#include "Components/TextBlock.h"
#include "Components/Slider.h"
#include "Components/CheckBox.h"
#include "Components/ComboBoxString.h"
#include "Kismet/GameplayStatics.h"

void UMDControlWidget::NativeConstruct()
{
    Super::NativeConstruct();
    
    // Find or spawn MD simulation actor
    FindMDSimulation();
    
    // Bind button events
    if (Button_Start)
        Button_Start->OnClicked.AddDynamic(this, &UMDControlWidget::OnStartClicked);
    if (Button_Stop)
        Button_Stop->OnClicked.AddDynamic(this, &UMDControlWidget::OnStopClicked);
    if (Button_Reset)
        Button_Reset->OnClicked.AddDynamic(this, &UMDControlWidget::OnResetClicked);
    if (Button_Initialize)
        Button_Initialize->OnClicked.AddDynamic(this, &UMDControlWidget::OnInitializeClicked);
    
    // Bind slider events
    if (Slider_TimeStep)
    {
        Slider_TimeStep->SetMinValue(0.0001f);
        Slider_TimeStep->SetMaxValue(0.01f);
        Slider_TimeStep->SetValue(0.001f);
        Slider_TimeStep->OnValueChanged.AddDynamic(this, &UMDControlWidget::OnTimeStepChanged);
    }
    
    if (Slider_Temperature)
    {
        Slider_Temperature->SetMinValue(0.0f);
        Slider_Temperature->SetMaxValue(1000.0f);
        Slider_Temperature->SetValue(300.0f);
        Slider_Temperature->OnValueChanged.AddDynamic(this, &UMDControlWidget::OnTemperatureChanged);
    }
    
    if (Slider_Damping)
    {
        Slider_Damping->SetMinValue(0.0f);
        Slider_Damping->SetMaxValue(1.0f);
        Slider_Damping->SetValue(0.95f);
        Slider_Damping->OnValueChanged.AddDynamic(this, &UMDControlWidget::OnDampingChanged);
    }
    
    // Bind checkbox events
    if (CheckBox_LennardJones)
    {
        CheckBox_LennardJones->SetIsChecked(true);
        CheckBox_LennardJones->OnCheckStateChanged.AddDynamic(this, &UMDControlWidget::OnLennardJonesChanged);
    }
    
    if (CheckBox_Electrostatics)
    {
        CheckBox_Electrostatics->SetIsChecked(false);
        CheckBox_Electrostatics->OnCheckStateChanged.AddDynamic(this, &UMDControlWidget::OnElectrostaticsChanged);
    }
    
    if (CheckBox_Constraints)
    {
        CheckBox_Constraints->SetIsChecked(true);
        CheckBox_Constraints->OnCheckStateChanged.AddDynamic(this, &UMDControlWidget::OnConstraintsChanged);
    }
    
    // Setup integrator combo box
    if (ComboBox_Integrator)
    {
        ComboBox_Integrator->ClearOptions();
        ComboBox_Integrator->AddOption(TEXT("Euler"));
        ComboBox_Integrator->AddOption(TEXT("Verlet"));
        ComboBox_Integrator->AddOption(TEXT("Leap-Frog"));
        ComboBox_Integrator->SetSelectedOption(TEXT("Verlet"));
        ComboBox_Integrator->OnSelectionChanged.AddDynamic(this, &UMDControlWidget::OnIntegratorChanged);
    }
    
    UpdateInfoDisplays();
}

void UMDControlWidget::NativeTick(const FGeometry& MyGeometry, float InDeltaTime)
{
    Super::NativeTick(MyGeometry, InDeltaTime);
    
    // Update info displays every frame
    UpdateInfoDisplays();
}

void UMDControlWidget::FindMDSimulation()
{
    // Try to find existing MD simulation
    TArray<AActor*> FoundActors;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), AMolecularDynamics::StaticClass(), FoundActors);
    
    if (FoundActors.Num() > 0)
    {
        MDSimulation = Cast<AMolecularDynamics>(FoundActors[0]);
        UE_LOG(LogTemp, Log, TEXT("MD Widget: Found existing MD simulation"));
    }
    else
    {
        // Spawn new MD simulation actor
        FActorSpawnParameters SpawnParams;
        SpawnParams.Name = FName(TEXT("MolecularDynamicsSimulation"));
        MDSimulation = GetWorld()->SpawnActor<AMolecularDynamics>(
            AMolecularDynamics::StaticClass(),
            FVector::ZeroVector,
            FRotator::ZeroRotator,
            SpawnParams
        );
        UE_LOG(LogTemp, Log, TEXT("MD Widget: Spawned new MD simulation"));
    }
}

void UMDControlWidget::OnStartClicked()
{
    if (MDSimulation)
    {
        MDSimulation->StartSimulation();
        UE_LOG(LogTemp, Log, TEXT("MD Widget: Started simulation"));
    }
}

void UMDControlWidget::OnStopClicked()
{
    if (MDSimulation)
    {
        MDSimulation->StopSimulation();
        UE_LOG(LogTemp, Log, TEXT("MD Widget: Stopped simulation"));
    }
}

void UMDControlWidget::OnResetClicked()
{
    if (MDSimulation)
    {
        MDSimulation->ResetSimulation();
        UE_LOG(LogTemp, Log, TEXT("MD Widget: Reset simulation"));
    }
}

void UMDControlWidget::OnInitializeClicked()
{
    if (!MDSimulation)
        return;
    
    // Find PDB Viewer
    TArray<AActor*> FoundActors;
    UGameplayStatics::GetAllActorsOfClass(GetWorld(), APDBViewer::StaticClass(), FoundActors);
    
    if (FoundActors.Num() > 0)
    {
        APDBViewer* Viewer = Cast<APDBViewer>(FoundActors[0]);
        if (Viewer)
        {
            // Count visible ligands before initialization
            int32 VisibleCount = 0;
            for (const auto& Pair : Viewer->LigandMap)
            {
                if (Pair.Value && Pair.Value->bIsVisible)
                    VisibleCount++;
            }
            
            if (VisibleCount == 0)
            {
                UE_LOG(LogTemp, Warning, TEXT("MD Widget: No visible molecules. Please make at least one molecule visible."));
            }
            else
            {
                MDSimulation->InitializeFromViewer(Viewer);
                UE_LOG(LogTemp, Log, TEXT("MD Widget: Initialized from %d visible molecule(s)"), VisibleCount);
            }
        }
    }
    else
    {
        UE_LOG(LogTemp, Warning, TEXT("MD Widget: No PDB Viewer found"));
    }
}

void UMDControlWidget::OnTimeStepChanged(float Value)
{
    if (MDSimulation)
    {
        MDSimulation->SetTimeStep(Value);
    }
    
    if (Text_TimeStep)
    {
        Text_TimeStep->SetText(FText::FromString(FString::Printf(TEXT("%.4f"), Value)));
    }
}

void UMDControlWidget::OnTemperatureChanged(float Value)
{
    if (MDSimulation)
    {
        MDSimulation->SetTemperature(Value);
    }
    
    if (Text_TargetTemp)
    {
        Text_TargetTemp->SetText(FText::FromString(FString::Printf(TEXT("%.1f K"), Value)));
    }
}

void UMDControlWidget::OnDampingChanged(float Value)
{
    if (MDSimulation)
    {
        MDSimulation->SetDamping(Value);
    }
    
    if (Text_DampingValue)
    {
        Text_DampingValue->SetText(FText::FromString(FString::Printf(TEXT("%.2f"), Value)));
    }
}

void UMDControlWidget::OnLennardJonesChanged(bool bIsChecked)
{
    // This would require adding a setter in AMolecularDynamics
    UE_LOG(LogTemp, Log, TEXT("MD Widget: Lennard-Jones %s"), 
        bIsChecked ? TEXT("enabled") : TEXT("disabled"));
}

void UMDControlWidget::OnElectrostaticsChanged(bool bIsChecked)
{
    UE_LOG(LogTemp, Log, TEXT("MD Widget: Electrostatics %s"), 
        bIsChecked ? TEXT("enabled") : TEXT("disabled"));
}

void UMDControlWidget::OnConstraintsChanged(bool bIsChecked)
{
    UE_LOG(LogTemp, Log, TEXT("MD Widget: Constraints %s"), 
        bIsChecked ? TEXT("enabled") : TEXT("disabled"));
}

void UMDControlWidget::OnIntegratorChanged(FString SelectedItem, ESelectInfo::Type SelectionType)
{
    if (!MDSimulation)
        return;
    
    EMDIntegrator NewIntegrator = EMDIntegrator::Verlet;
    
    if (SelectedItem == TEXT("Euler"))
        NewIntegrator = EMDIntegrator::Euler;
    else if (SelectedItem == TEXT("Verlet"))
        NewIntegrator = EMDIntegrator::Verlet;
    else if (SelectedItem == TEXT("Leap-Frog"))
        NewIntegrator = EMDIntegrator::LeapFrog;
    
    MDSimulation->SetIntegrator(NewIntegrator);
    UE_LOG(LogTemp, Log, TEXT("MD Widget: Set integrator to %s"), *SelectedItem);
}

void UMDControlWidget::UpdateInfoDisplays()
{
    if (!MDSimulation)
        return;
    
    // Update status
    if (Text_Status)
    {
        FString Status = MDSimulation->IsSimulating() ? TEXT("Running") : TEXT("Stopped");
        Text_Status->SetText(FText::FromString(Status));
        
        // Color coding
        FLinearColor StatusColor = MDSimulation->IsSimulating() ? 
            FLinearColor::Green : FLinearColor::Gray;
        Text_Status->SetColorAndOpacity(FSlateColor(StatusColor));
    }
    
    // Update atom count
    if (Text_AtomCount)
    {
        int32 Count = MDSimulation->GetAtomCount();
        Text_AtomCount->SetText(FText::FromString(FString::Printf(TEXT("%d atoms"), Count)));
    }
    
    // Update energy
    if (Text_Energy)
    {
        float TotalE = MDSimulation->GetCurrentEnergy();
        float KineticE = MDSimulation->GetKineticEnergy();
        float PotentialE = MDSimulation->GetPotentialEnergy();
        
        FString EnergyText = FString::Printf(
            TEXT("Total: %.2f\nKinetic: %.2f\nPotential: %.2f"),
            TotalE, KineticE, PotentialE
        );
        Text_Energy->SetText(FText::FromString(EnergyText));
    }
    
    // Update temperature (calculated from kinetic energy)
    if (Text_Temperature)
    {
        float KE = MDSimulation->GetKineticEnergy();
        int32 N = MDSimulation->GetAtomCount();
        
        // Simplified temperature from KE
        float CurrentTemp = N > 0 ? KE / (N * 1.5f) : 0.0f;
        
        Text_Temperature->SetText(FText::FromString(
            FString::Printf(TEXT("%.1f K"), CurrentTemp)
        ));
    }
}