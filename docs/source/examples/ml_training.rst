ML Model Training
===============

This guide shows how to train machine learning models for various prediction tasks.

Dataset Preparation
----------------

Preparing training data:

.. code-block:: python

    from binding_data_processor.pipeline.ml import DatasetBuilder
    from binding_data_processor.models import CompoundData
    from typing import List, Tuple

    def prepare_dataset(
        compounds: List[CompoundData],
        split_ratio: float = 0.8,
    ) -> Tuple[List[CompoundData], List[CompoundData]]:
        """Prepare training and validation datasets."""
        
        # Configure builder
        builder = DatasetBuilder(
            feature_types=["fingerprints", "descriptors"],
            label_types=["binding", "activity"],
            split_ratio=split_ratio,
        )
        
        # Build datasets
        train_data, valid_data = builder.build_datasets(compounds)
        
        return train_data, valid_data

Binding Predictor
--------------

Training a binding affinity predictor:

.. code-block:: python

    import torch
    import torch.nn as nn
    from binding_data_processor.pipeline.ml import (
        BindingPredictor,
        ModelTrainer,
        TrainingConfig,
    )

    # Define model architecture
    class BindingModel(nn.Module):
        def __init__(self, input_size: int, num_receptors: int):
            super().__init__()
            self.layers = nn.Sequential(
                nn.Linear(input_size, 512),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(512, 256),
                nn.ReLU(),
                nn.Linear(256, num_receptors),
                nn.Sigmoid(),
            )

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            return self.layers(x)

    # Configure training
    config = TrainingConfig(
        batch_size=32,
        learning_rate=0.001,
        num_epochs=100,
        early_stopping=True,
        patience=10,
    )

    # Create trainer
    trainer = ModelTrainer(
        model=BindingModel(input_size=1024, num_receptors=5),
        config=config,
    )

    # Train model
    model = trainer.train(
        train_data=train_data,
        valid_data=valid_data,
        model_path="models/binding/model.pt",
    )

Activity Predictor
---------------

Training an activity predictor:

.. code-block:: python

    from binding_data_processor.pipeline.ml import ActivityPredictor

    class ActivityModel(nn.Module):
        def __init__(self, input_size: int, num_effects: int):
            super().__init__()
            self.layers = nn.Sequential(
                nn.Linear(input_size, 512),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(512, 256),
                nn.ReLU(),
                nn.Linear(256, num_effects),
                nn.Sigmoid(),
            )

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            return self.layers(x)

    # Configure training
    config = TrainingConfig(
        batch_size=32,
        learning_rate=0.001,
        num_epochs=100,
        early_stopping=True,
        patience=10,
    )

    # Create trainer
    trainer = ModelTrainer(
        model=ActivityModel(input_size=1024, num_effects=10),
        config=config,
    )

    # Train model
    model = trainer.train(
        train_data=train_data,
        valid_data=valid_data,
        model_path="models/activity/model.pt",
    )

Safety Predictor
-------------

Training a safety predictor:

.. code-block:: python

    from binding_data_processor.pipeline.ml import SafetyPredictor

    class SafetyModel(nn.Module):
        def __init__(self, input_size: int, num_risks: int):
            super().__init__()
            self.layers = nn.Sequential(
                nn.Linear(input_size, 512),
                nn.ReLU(),
                nn.Dropout(0.5),
                nn.Linear(512, 256),
                nn.ReLU(),
                nn.Linear(256, num_risks),
                nn.Sigmoid(),
            )

        def forward(self, x: torch.Tensor) -> torch.Tensor:
            return self.layers(x)

    # Configure training
    config = TrainingConfig(
        batch_size=32,
        learning_rate=0.001,
        num_epochs=100,
        early_stopping=True,
        patience=10,
    )

    # Create trainer
    trainer = ModelTrainer(
        model=SafetyModel(input_size=1024, num_risks=8),
        config=config,
    )

    # Train model
    model = trainer.train(
        train_data=train_data,
        valid_data=valid_data,
        model_path="models/safety/model.pt",
    )

Model Ensembles
------------

Creating model ensembles:

.. code-block:: python

    from binding_data_processor.pipeline.ml import EnsembleTrainer

    # Configure ensemble
    config = EnsembleConfig(
        num_models=5,
        model_type="binding",
        aggregation="weighted_average",
        weights=[0.3, 0.2, 0.2, 0.15, 0.15],
    )

    # Create trainer
    trainer = EnsembleTrainer(
        base_model=BindingModel(input_size=1024, num_receptors=5),
        config=config,
    )

    # Train ensemble
    ensemble = trainer.train(
        train_data=train_data,
        valid_data=valid_data,
        model_dir="models/binding/ensemble/",
    )

Transfer Learning
--------------

Fine-tuning pre-trained models:

.. code-block:: python

    from binding_data_processor.pipeline.ml import TransferLearner

    # Configure transfer learning
    config = TransferConfig(
        base_model_path="models/binding/pretrained.pt",
        freeze_layers=True,
        learning_rate=0.0001,
        num_epochs=50,
    )

    # Create learner
    learner = TransferLearner(
        model=BindingModel(input_size=1024, num_receptors=5),
        config=config,
    )

    # Fine-tune model
    model = learner.train(
        train_data=train_data,
        valid_data=valid_data,
        model_path="models/binding/finetuned.pt",
    )

Model Evaluation
-------------

Evaluating model performance:

.. code-block:: python

    from binding_data_processor.pipeline.ml import ModelEvaluator

    # Configure evaluator
    evaluator = ModelEvaluator(
        metrics=["accuracy", "precision", "recall", "f1"],
        num_folds=5,
    )

    # Evaluate model
    results = evaluator.evaluate(
        model=model,
        test_data=test_data,
    )

    # Print results
    for metric, value in results.items():
        print(f"{metric}: {value:.3f}")

    # Generate confusion matrix
    evaluator.plot_confusion_matrix(
        model=model,
        test_data=test_data,
        output_file="confusion_matrix.png",
    )

Hyperparameter Tuning
------------------

Optimizing model hyperparameters:

.. code-block:: python

    from binding_data_processor.pipeline.ml import HyperparameterTuner

    # Define parameter space
    param_space = {
        "learning_rate": [0.1, 0.01, 0.001],
        "batch_size": [16, 32, 64],
        "hidden_size": [256, 512, 1024],
        "dropout": [0.3, 0.5, 0.7],
    }

    # Configure tuner
    tuner = HyperparameterTuner(
        model_class=BindingModel,
        param_space=param_space,
        num_trials=50,
        cv_folds=5,
    )

    # Find best parameters
    best_params = tuner.optimize(
        train_data=train_data,
        valid_data=valid_data,
    )

    print(f"Best parameters: {best_params}")

Model Deployment
-------------

Saving and loading models:

.. code-block:: python

    # Save model
    trainer.save_model(
        model=model,
        path="models/binding/production.pt",
        include_config=True,
    )

    # Load model
    loaded_model = trainer.load_model(
        path="models/binding/production.pt",
    )

    # Create predictor
    predictor = BindingPredictor(
        model=loaded_model,
        config=predictor_config,
    )

    # Make predictions
    for compound in compounds:
        prediction = predictor.predict(compound)
        print(f"Prediction: {prediction.value}")
        print(f"Confidence: {prediction.confidence}")

Model Monitoring
-------------

Monitoring model performance:

.. code-block:: python

    from binding_data_processor.pipeline.ml import ModelMonitor

    # Configure monitor
    monitor = ModelMonitor(
        model=model,
        metrics=["accuracy", "drift"],
        alert_threshold=0.1,
    )

    # Monitor predictions
    for batch in data_stream:
        metrics = monitor.track_predictions(batch)
        if monitor.should_retrain():
            print("Model retraining needed")

    # Generate monitoring report
    report = monitor.generate_report()
    report.save("monitoring_report.html")
