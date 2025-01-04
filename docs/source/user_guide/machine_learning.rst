Machine Learning Guide
===================

This guide covers the machine learning capabilities of ChemData in detail.

ML Pipeline Overview
-----------------

The ML pipeline includes several predictors:

1. Binding Affinity Prediction
2. BBB Permeability Prediction
3. Activity Classification
4. Toxicity Prediction
5. Abuse Potential Prediction

Basic Usage
---------

Using pre-trained models:

.. code-block:: python

    from binding_data_processor.pipeline import ProcessingPipeline
    from binding_data_processor.pipeline.config import ProcessingConfig

    # Create pipeline with ML enabled
    pipeline = ProcessingPipeline(
        config=ProcessingConfig(
            use_ml_predictions=True,
            model_dir="models/",
        )
    )

    # Process compounds
    compounds = pipeline.process_compounds(
        input_file="compounds.tsv",
        output_dir="results/",
    )

Binding Prediction
---------------

Predicting binding affinities:

.. code-block:: python

    from binding_data_processor.processors.psychopharm.predictors import (
        BindingPredictor,
        BindingConfig,
    )

    # Configure predictor
    predictor = BindingPredictor(
        config=BindingConfig(
            target_receptors=["5-HT2A", "NMDA"],
            confidence_threshold=0.8,
        )
    )

    # Make predictions
    for compound in compounds:
        predictions = predictor.predict(compound)
        print(f"Compound: {compound.name}")
        for receptor, (affinity, confidence) in predictions.items():
            print(f"- {receptor}: {affinity:.2f} (conf: {confidence:.2f})")

BBB Prediction
-----------

Predicting blood-brain barrier permeability:

.. code-block:: python

    from binding_data_processor.processors.psychopharm.predictors.bbb import (
        BBBPredictorWebEnriched,
        BBBConfig,
    )

    # Configure predictor
    predictor = BBBPredictorWebEnriched(
        config=BBBConfig(
            use_transporter_data=True,
            use_web_data=True,
        )
    )

    # Make predictions
    for compound in compounds:
        result = predictor.predict(compound)
        print(f"Compound: {compound.name}")
        print(f"BBB Class: {result.value}")
        print(f"Confidence: {result.confidence:.2f}")
        print("Supporting Data:")
        for key, value in result.supporting_data.items():
            print(f"- {key}: {value}")

Activity Prediction
----------------

Predicting psychoactive activity:

.. code-block:: python

    from binding_data_processor.processors.psychopharm.predictors import (
        ActivityPredictor,
        ActivityConfig,
    )

    # Configure predictor
    predictor = ActivityPredictor(
        config=ActivityConfig(
            effect_types=["stimulant", "psychedelic", "depressant"],
            min_confidence=0.7,
        )
    )

    # Make predictions
    for compound in compounds:
        predictions = predictor.predict(compound)
        print(f"Compound: {compound.name}")
        for effect, (score, confidence) in predictions.items():
            print(f"- {effect}: {score:.2f} (conf: {confidence:.2f})")

Toxicity Prediction
----------------

Predicting toxicity:

.. code-block:: python

    from binding_data_processor.processors.psychopharm.predictors import (
        ToxicityPredictor,
        ToxicityConfig,
    )

    # Configure predictor
    predictor = ToxicityPredictor(
        config=ToxicityConfig(
            endpoints=["acute", "chronic", "cardio", "hepato"],
            threshold=0.5,
        )
    )

    # Make predictions
    for compound in compounds:
        predictions = predictor.predict(compound)
        print(f"Compound: {compound.name}")
        for endpoint, (risk, confidence) in predictions.items():
            print(f"- {endpoint}: {risk:.2f} (conf: {confidence:.2f})")

Abuse Prediction
-------------

Predicting abuse potential:

.. code-block:: python

    from binding_data_processor.processors.psychopharm.predictors import (
        AbusePredictor,
        AbuseConfig,
    )

    # Configure predictor
    predictor = AbusePredictor(
        config=AbuseConfig(
            use_binding_data=True,
            use_community_data=True,
        )
    )

    # Make predictions
    for compound in compounds:
        result = predictor.predict(compound)
        print(f"Compound: {compound.name}")
        print(f"Abuse Potential: {result.score:.2f}")
        print(f"Confidence: {result.confidence:.2f}")
        print("Risk Factors:")
        for factor, score in result.risk_factors.items():
            print(f"- {factor}: {score:.2f}")

Model Training
------------

Training Custom Models
~~~~~~~~~~~~~~~~~~

Training a new model:

.. code-block:: python

    from binding_data_processor.pipeline.ml import ModelTrainer

    # Configure trainer
    trainer = ModelTrainer(
        model_type="binding",
        target_receptors=["5-HT2A"],
        hyperparameters={
            "learning_rate": 0.001,
            "batch_size": 32,
            "epochs": 100,
        },
    )

    # Train model
    model = trainer.train(
        train_data="train.tsv",
        valid_data="valid.tsv",
        output_dir="models/binding/",
    )

Transfer Learning
~~~~~~~~~~~~~~

Fine-tuning existing models:

.. code-block:: python

    from binding_data_processor.pipeline.ml import TransferLearner

    # Configure learner
    learner = TransferLearner(
        base_model="models/binding/5ht2a.pt",
        freeze_layers=True,
        learning_rate=0.0001,
    )

    # Fine-tune model
    model = learner.train(
        train_data="new_data.tsv",
        valid_data="new_valid.tsv",
        output_dir="models/binding/custom/",
    )

Ensemble Methods
-------------

Using model ensembles:

.. code-block:: python

    from binding_data_processor.pipeline.ml import EnsemblePredictor

    # Configure ensemble
    ensemble = EnsemblePredictor(
        models=[
            "models/binding/model1.pt",
            "models/binding/model2.pt",
            "models/binding/model3.pt",
        ],
        weights=[0.4, 0.3, 0.3],
    )

    # Make predictions
    for compound in compounds:
        prediction = ensemble.predict(compound)
        print(f"Ensemble Prediction: {prediction.value}")
        print(f"Confidence: {prediction.confidence:.2f}")

Uncertainty Estimation
------------------

Estimating prediction uncertainty:

.. code-block:: python

    from binding_data_processor.pipeline.ml import UncertaintyEstimator

    # Configure estimator
    estimator = UncertaintyEstimator(
        method="dropout",
        num_samples=100,
    )

    # Get uncertainty estimates
    for compound in compounds:
        prediction = predictor.predict(compound)
        uncertainty = estimator.estimate(prediction)
        print(f"Prediction: {prediction.value}")
        print(f"Uncertainty: {uncertainty:.2f}")

Model Validation
-------------

Validating model performance:

.. code-block:: python

    from binding_data_processor.pipeline.ml import ModelValidator

    # Configure validator
    validator = ModelValidator(
        metrics=["accuracy", "precision", "recall", "f1"],
        cross_validation=True,
        num_folds=5,
    )

    # Validate model
    results = validator.validate(
        model=predictor,
        test_data="test.tsv",
    )

    # Print results
    for metric, value in results.items():
        print(f"{metric}: {value:.3f}")

Advanced Usage
------------

Custom Features
~~~~~~~~~~~~

Creating custom molecular features:

.. code-block:: python

    from binding_data_processor.pipeline.ml import FeatureExtractor

    class CustomExtractor(FeatureExtractor):
        def extract_features(self, compound):
            # Custom feature extraction
            features = {
                "custom_fp": self.calc_custom_fingerprint(compound),
                "custom_desc": self.calc_custom_descriptors(compound),
            }
            return features

    # Use custom extractor
    extractor = CustomExtractor()
    features = extractor.extract_features(compound)

Custom Models
~~~~~~~~~~~

Creating custom ML models:

.. code-block:: python

    from binding_data_processor.pipeline.ml import BaseModel
    import torch.nn as nn

    class CustomModel(BaseModel):
        def __init__(self):
            super().__init__()
            self.layers = nn.Sequential(
                nn.Linear(1024, 512),
                nn.ReLU(),
                nn.Linear(512, 1),
            )

        def forward(self, x):
            return self.layers(x)

    # Use custom model
    model = CustomModel()
    trainer.train(model, train_data)
