#!/bin/bash
# Script to manage machine learning models, training, and evaluation

# Exit on error
set -e

# Default values
DATA_DIR="data/models"
MODEL_DIR="models"
CACHE_DIR="cache/models"
LOG_DIR="logs"
MODEL_TYPES="activity,toxicity,abuse,psychoactive,nootropic,bbb,binding,affinity,ensemble"
MODEL_STAGES="base,tuned,ensemble"
MODEL_MODES="development,production,testing"
DATASET_TYPES="bindingdb,chembl,pubchem,swiss,community,literature"
BATCH_SIZES="32,64,128,256"
EPOCHS="100,200,500,1000"
LEARNING_RATES="0.1,0.01,0.001,0.0001"
TRAIN=false
TUNE=false
EVALUATE=false
PREDICT=false
OPTIMIZE=false
EXPORT=false
CLEAN=false
BACKUP=false
RESTORE=false
GPU=false
FORCE=false

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --model-dir)
            MODEL_DIR="$2"
            shift 2
            ;;
        --cache-dir)
            CACHE_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;
        --model-types)
            MODEL_TYPES="$2"
            shift 2
            ;;
        --model-stages)
            MODEL_STAGES="$2"
            shift 2
            ;;
        --model-modes)
            MODEL_MODES="$2"
            shift 2
            ;;
        --dataset-types)
            DATASET_TYPES="$2"
            shift 2
            ;;
        --batch-sizes)
            BATCH_SIZES="$2"
            shift 2
            ;;
        --epochs)
            EPOCHS="$2"
            shift 2
            ;;
        --learning-rates)
            LEARNING_RATES="$2"
            shift 2
            ;;
        --train)
            TRAIN=true
            shift
            ;;
        --tune)
            TUNE=true
            shift
            ;;
        --evaluate)
            EVALUATE=true
            shift
            ;;
        --predict)
            PREDICT=true
            shift
            ;;
        --optimize)
            OPTIMIZE=true
            shift
            ;;
        --export)
            EXPORT=true
            shift
            ;;
        --clean)
            CLEAN=true
            shift
            ;;
        --backup)
            BACKUP=true
            shift
            ;;
        --restore)
            RESTORE=true
            shift
            ;;
        --gpu)
            GPU=true
            shift
            ;;
        --force)
            FORCE=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to create directory structure
create_dirs() {
    echo "Creating directory structure..."
    
    # Model directories
    for type in ${MODEL_TYPES//,/ }; do
        for stage in ${MODEL_STAGES//,/ }; do
            for mode in ${MODEL_MODES//,/ }; do
                mkdir -p "$MODEL_DIR/$type/$stage/$mode"/{checkpoints,configs,metrics,reports}
            done
        done
    done
    
    # Data directories
    for type in ${DATASET_TYPES//,/ }; do
        mkdir -p "$DATA_DIR/$type"/{raw,processed,features,splits}
    done
    
    # Cache directory
    mkdir -p "$CACHE_DIR"
    
    # Log directory
    mkdir -p "$LOG_DIR"
}

# Function to train models
train_models() {
    echo "Training models..."
    
    for mode in ${MODEL_MODES//,/ }; do
        echo "Training models in $mode mode..."
        
        for type in ${MODEL_TYPES//,/ }; do
            echo "Training $type models..."
            
            for stage in ${MODEL_STAGES//,/ }; do
                echo "Training $stage model..."
                
                # Train on each dataset type
                for dataset in ${DATASET_TYPES//,/ }; do
                    echo "Training on $dataset dataset..."
                    
                    # Try different hyperparameters
                    for batch_size in ${BATCH_SIZES//,/ }; do
                        for epochs in ${EPOCHS//,/ }; do
                            for lr in ${LEARNING_RATES//,/ }; do
                                echo "Training with batch_size=$batch_size, epochs=$epochs, lr=$lr..."
                                
                                # Build training command
                                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli train-model"
                                CMD="$CMD --data-dir $DATA_DIR/$dataset"
                                CMD="$CMD --model-dir $MODEL_DIR/$type/$stage/$mode"
                                CMD="$CMD --cache-dir $CACHE_DIR"
                                CMD="$CMD --type $type"
                                CMD="$CMD --stage $stage"
                                CMD="$CMD --mode $mode"
                                CMD="$CMD --dataset $dataset"
                                CMD="$CMD --batch-size $batch_size"
                                CMD="$CMD --epochs $epochs"
                                CMD="$CMD --learning-rate $lr"
                                CMD="$CMD --log-dir $LOG_DIR"
                                
                                if [ "$GPU" = true ]; then
                                    CMD="$CMD --gpu"
                                fi
                                
                                # Run training
                                echo "Running: $CMD"
                                $CMD || true  # Continue even if training fails
                                
                                # Generate training report
                                echo "Generating training report..."
                                python -m binding_data_processor.processors.psychopharm.predictors.cli generate-training-report \
                                    --model-dir "$MODEL_DIR/$type/$stage/$mode" \
                                    --dataset "$dataset" \
                                    --batch-size "$batch_size" \
                                    --epochs "$epochs" \
                                    --learning-rate "$lr" \
                                    --output "$MODEL_DIR/$type/$stage/$mode/reports/training_${dataset}_${batch_size}_${epochs}_${lr}.html"
                            done
                        done
                    done
                done
            done
        done
    done
}

# Function to tune models
tune_models() {
    echo "Tuning models..."
    
    for mode in ${MODEL_MODES//,/ }; do
        echo "Tuning models in $mode mode..."
        
        for type in ${MODEL_TYPES//,/ }; do
            echo "Tuning $type models..."
            
            # Build tuning command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli tune-model"
            CMD="$CMD --data-dir $DATA_DIR"
            CMD="$CMD --model-dir $MODEL_DIR/$type/base/$mode"
            CMD="$CMD --cache-dir $CACHE_DIR"
            CMD="$CMD --type $type"
            CMD="$CMD --mode $mode"
            CMD="$CMD --batch-sizes $BATCH_SIZES"
            CMD="$CMD --epochs $EPOCHS"
            CMD="$CMD --learning-rates $LEARNING_RATES"
            CMD="$CMD --log-dir $LOG_DIR"
            
            if [ "$GPU" = true ]; then
                CMD="$CMD --gpu"
            fi
            
            # Run tuning
            echo "Running: $CMD"
            $CMD || true  # Continue even if tuning fails
            
            # Generate tuning report
            echo "Generating tuning report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-tuning-report \
                --model-dir "$MODEL_DIR/$type/tuned/$mode" \
                --output "$MODEL_DIR/$type/tuned/$mode/reports/tuning.html"
        done
    done
}

# Function to evaluate models
evaluate_models() {
    echo "Evaluating models..."
    
    for mode in ${MODEL_MODES//,/ }; do
        echo "Evaluating models in $mode mode..."
        
        for type in ${MODEL_TYPES//,/ }; do
            echo "Evaluating $type models..."
            
            for stage in ${MODEL_STAGES//,/ }; do
                echo "Evaluating $stage model..."
                
                # Evaluate on each dataset type
                for dataset in ${DATASET_TYPES//,/ }; do
                    echo "Evaluating on $dataset dataset..."
                    
                    # Build evaluation command
                    CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli evaluate-model"
                    CMD="$CMD --data-dir $DATA_DIR/$dataset"
                    CMD="$CMD --model-dir $MODEL_DIR/$type/$stage/$mode"
                    CMD="$CMD --cache-dir $CACHE_DIR"
                    CMD="$CMD --type $type"
                    CMD="$CMD --stage $stage"
                    CMD="$CMD --mode $mode"
                    CMD="$CMD --dataset $dataset"
                    CMD="$CMD --log-dir $LOG_DIR"
                    
                    if [ "$GPU" = true ]; then
                        CMD="$CMD --gpu"
                    fi
                    
                    # Run evaluation
                    echo "Running: $CMD"
                    $CMD || true  # Continue even if evaluation fails
                    
                    # Generate evaluation report
                    echo "Generating evaluation report..."
                    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-evaluation-report \
                        --model-dir "$MODEL_DIR/$type/$stage/$mode" \
                        --dataset "$dataset" \
                        --output "$MODEL_DIR/$type/$stage/$mode/reports/evaluation_${dataset}.html"
                done
            done
        done
    done
}

# Function to make predictions
make_predictions() {
    echo "Making predictions..."
    
    for mode in ${MODEL_MODES//,/ }; do
        echo "Making predictions in $mode mode..."
        
        for type in ${MODEL_TYPES//,/ }; do
            echo "Making $type predictions..."
            
            for stage in ${MODEL_STAGES//,/ }; do
                echo "Making predictions with $stage model..."
                
                # Predict on each dataset type
                for dataset in ${DATASET_TYPES//,/ }; do
                    echo "Predicting on $dataset dataset..."
                    
                    # Build prediction command
                    CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli predict"
                    CMD="$CMD --data-dir $DATA_DIR/$dataset"
                    CMD="$CMD --model-dir $MODEL_DIR/$type/$stage/$mode"
                    CMD="$CMD --cache-dir $CACHE_DIR"
                    CMD="$CMD --type $type"
                    CMD="$CMD --stage $stage"
                    CMD="$CMD --mode $mode"
                    CMD="$CMD --dataset $dataset"
                    CMD="$CMD --log-dir $LOG_DIR"
                    
                    if [ "$GPU" = true ]; then
                        CMD="$CMD --gpu"
                    fi
                    
                    # Run prediction
                    echo "Running: $CMD"
                    $CMD || true  # Continue even if prediction fails
                    
                    # Generate prediction report
                    echo "Generating prediction report..."
                    python -m binding_data_processor.processors.psychopharm.predictors.cli generate-prediction-report \
                        --model-dir "$MODEL_DIR/$type/$stage/$mode" \
                        --dataset "$dataset" \
                        --output "$MODEL_DIR/$type/$stage/$mode/reports/prediction_${dataset}.html"
                done
            done
        done
    done
}

# Function to optimize models
optimize_models() {
    echo "Optimizing models..."
    
    for mode in ${MODEL_MODES//,/ }; do
        echo "Optimizing models in $mode mode..."
        
        for type in ${MODEL_TYPES//,/ }; do
            echo "Optimizing $type models..."
            
            # Optimize on each dataset type
            for dataset in ${DATASET_TYPES//,/ }; do
                echo "Optimizing on $dataset dataset..."
                
                # Build optimization command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli optimize-model"
                CMD="$CMD --data-dir $DATA_DIR/$dataset"
                CMD="$CMD --model-dir $MODEL_DIR/$type/base/$mode"
                CMD="$CMD --cache-dir $CACHE_DIR"
                CMD="$CMD --type $type"
                CMD="$CMD --mode $mode"
                CMD="$CMD --dataset $dataset"
                CMD="$CMD --log-dir $LOG_DIR"
                
                if [ "$GPU" = true ]; then
                    CMD="$CMD --gpu"
                fi
                
                # Run optimization
                echo "Running: $CMD"
                $CMD || true  # Continue even if optimization fails
                
                # Generate optimization report
                echo "Generating optimization report..."
                python -m binding_data_processor.processors.psychopharm.predictors.cli generate-optimization-report \
                    --model-dir "$MODEL_DIR/$type/tuned/$mode" \
                    --dataset "$dataset" \
                    --output "$MODEL_DIR/$type/tuned/$mode/reports/optimization_${dataset}.html"
            done
        done
    done
}

# Function to export models
export_models() {
    echo "Exporting models..."
    
    for mode in ${MODEL_MODES//,/ }; do
        echo "Exporting models in $mode mode..."
        
        for type in ${MODEL_TYPES//,/ }; do
            echo "Exporting $type models..."
            
            for stage in ${MODEL_STAGES//,/ }; do
                echo "Exporting $stage model..."
                
                # Build export command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli export-model"
                CMD="$CMD --model-dir $MODEL_DIR/$type/$stage/$mode"
                CMD="$CMD --type $type"
                CMD="$CMD --stage $stage"
                CMD="$CMD --mode $mode"
                CMD="$CMD --log-dir $LOG_DIR"
                
                # Run export
                echo "Running: $CMD"
                $CMD || true  # Continue even if export fails
            done
        done
    done
}

# Function to clean models
clean_models() {
    echo "Cleaning models..."
    
    # Clean model directories
    rm -rf "$MODEL_DIR"/*
    
    # Clean data directories
    rm -rf "$DATA_DIR"/*
    
    # Clean cache
    rm -rf "$CACHE_DIR"/*
}

# Function to backup models
backup_models() {
    echo "Backing up models..."
    
    # Create backup directory
    BACKUP_DIR="backups/models_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$MODEL_DIR" "$BACKUP_DIR/"
    cp -r "$DATA_DIR" "$BACKUP_DIR/"
    
    # Create encrypted backup if gpg is available
    if command_exists gpg; then
        echo "Encrypting backup..."
        tar -czf - "$BACKUP_DIR" | gpg --symmetric --output "$BACKUP_DIR.tar.gz.gpg"
        rm -rf "$BACKUP_DIR"
        echo "Encrypted backup saved to: $BACKUP_DIR.tar.gz.gpg"
    else
        # Create unencrypted backup
        tar -czf "$BACKUP_DIR.tar.gz" "$BACKUP_DIR"
        rm -rf "$BACKUP_DIR"
        echo "Backup saved to: $BACKUP_DIR.tar.gz"
    fi
}

# Function to restore models
restore_models() {
    echo "Restoring models..."
    
    # Find latest backup
    local LATEST_BACKUP=""
    if command_exists gpg; then
        LATEST_BACKUP=$(ls -t backups/models_*.tar.gz.gpg 2>/dev/null | head -n1)
        if [ -n "$LATEST_BACKUP" ]; then
            echo "Decrypting backup: $LATEST_BACKUP"
            gpg --decrypt "$LATEST_BACKUP" | tar -xz
            BACKUP_DIR="${LATEST_BACKUP%.tar.gz.gpg}"
        fi
    fi
    
    if [ -z "$LATEST_BACKUP" ]; then
        LATEST_BACKUP=$(ls -t backups/models_*.tar.gz 2>/dev/null | head -n1)
        if [ -n "$LATEST_BACKUP" ]; then
            echo "Extracting backup: $LATEST_BACKUP"
            tar -xzf "$LATEST_BACKUP"
            BACKUP_DIR="${LATEST_BACKUP%.tar.gz}"
        fi
    fi
    
    if [ -z "$BACKUP_DIR" ]; then
        echo "No backup found"
        exit 1
    fi
    
    # Restore directories
    if [ "$FORCE" = true ]; then
        rm -rf "$MODEL_DIR" "$DATA_DIR"
    fi
    
    cp -r "$BACKUP_DIR/models" ./
    cp -r "$BACKUP_DIR/data/models" ./data/
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Models restored from: $BACKUP_DIR"
}

# Function to show model statistics
show_stats() {
    echo "Model statistics:"
    echo
    
    echo "Model files:"
    for mode in ${MODEL_MODES//,/ }; do
        echo "$mode mode:"
        for type in ${MODEL_TYPES//,/ }; do
            echo "  $type models:"
            
            for stage in ${MODEL_STAGES//,/ }; do
                echo "    $stage model:"
                
                # Checkpoint files
                echo "      Checkpoints:"
                echo "        Files: $(find "$MODEL_DIR/$type/$stage/$mode/checkpoints" -type f | wc -l) files"
                echo "        Size: $(du -sh "$MODEL_DIR/$type/$stage/$mode/checkpoints" | cut -f1)"
                
                # Config files
                echo "      Configs:"
                echo "        Files: $(find "$MODEL_DIR/$type/$stage/$mode/configs" -type f | wc -l) files"
                echo "        Size: $(du -sh "$MODEL_DIR/$type/$stage/$mode/configs" | cut -f1)"
                
                # Metric files
                echo "      Metrics:"
                echo "        Files: $(find "$MODEL_DIR/$type/$stage/$mode/metrics" -type f | wc -l) files"
                echo "        Size: $(du -sh "$MODEL_DIR/$type/$stage/$mode/metrics" | cut -f1)"
                
                # Show dataset stats
                for dataset in ${DATASET_TYPES//,/ }; do
                    if [ -f "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json" ]; then
                        echo "      $dataset metrics:"
                        case "$type" in
                            activity)
                                echo "        Accuracy: $(jq '.accuracy' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        ROC AUC: $(jq '.roc_auc' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                ;;
                            toxicity)
                                echo "        Precision: $(jq '.precision' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        Recall: $(jq '.recall' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                ;;
                            abuse)
                                echo "        F1 Score: $(jq '.f1' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        MCC: $(jq '.mcc' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                ;;
                            psychoactive)
                                echo "        Accuracy: $(jq '.accuracy' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        Kappa: $(jq '.kappa' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                ;;
                            nootropic)
                                echo "        Accuracy: $(jq '.accuracy' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        F1 Score: $(jq '.f1' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                ;;
                            bbb)
                                echo "        Accuracy: $(jq '.accuracy' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        MCC: $(jq '.mcc' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                ;;
                            binding)
                                echo "        RMSE: $(jq '.rmse' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        R2: $(jq '.r2' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                ;;
                            affinity)
                                echo "        MAE: $(jq '.mae' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        Pearson: $(jq '.pearson' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                ;;
                            ensemble)
                                echo "        Accuracy: $(jq '.accuracy' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        F1 Score: $(jq '.f1' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                echo "        ROC AUC: $(jq '.roc_auc' "$MODEL_DIR/$type/$stage/$mode/metrics/${dataset}_metrics.json")"
                                ;;
                        esac
                    fi
                done
                
                # Show reports
                if [ -f "$MODEL_DIR/$type/$stage/$mode/reports/training.html" ]; then
                    echo "      Training report: $MODEL_DIR/$type/$stage/$mode/reports/training.html"
                fi
                if [ -f "$MODEL_DIR/$type/$stage/$mode/reports/tuning.html" ]; then
                    echo "      Tuning report: $MODEL_DIR/$type/$stage/$mode/reports/tuning.html"
                fi
                if [ -f "$MODEL_DIR/$type/$stage/$mode/reports/evaluation.html" ]; then
                    echo "      Evaluation report: $MODEL_DIR/$type/$stage/$mode/reports/evaluation.html"
                fi
                if [ -f "$MODEL_DIR/$type/$stage/$mode/reports/prediction.html" ]; then
                    echo "      Prediction report: $MODEL_DIR/$type/$stage/$mode/reports/prediction.html"
                fi
                if [ -f "$MODEL_DIR/$type/$stage/$mode/reports/optimization.html" ]; then
                    echo "      Optimization report: $MODEL_DIR/$type/$stage/$mode/reports/optimization.html"
                fi
                echo
            done
        done
    done
    
    echo "Data files:"
    for type in ${DATASET_TYPES//,/ }; do
        echo "$type dataset:"
        
        # Raw data
        echo "  Raw data:"
        echo "    Files: $(find "$DATA_DIR/$type/raw" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DATA_DIR/$type/raw" | cut -f1)"
        
        # Processed data
        echo "  Processed data:"
        echo "    Files: $(find "$DATA_DIR/$type/processed" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DATA_DIR/$type/processed" | cut -f1)"
        
        # Feature data
        echo "  Feature data:"
        echo "    Files: $(find "$DATA_DIR/$type/features" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DATA_DIR/$type/features" | cut -f1)"
        
        # Split data
        echo "  Split data:"
        echo "    Files: $(find "$DATA_DIR/$type/splits" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DATA_DIR/$type/splits" | cut -f1)"
        echo
    done
    
    echo "Cache:"
    echo "  Size: $(du -sh "$CACHE_DIR" | cut -f1)"
    echo "  Files: $(find "$CACHE_DIR" -type f | wc -l) files"
    echo
    
    echo "Logs:"
    echo "  Size: $(du -sh "$LOG_DIR" | cut -f1)"
    echo "  Files: $(find "$LOG_DIR" -type f | wc -l) files"
}

# Main process
echo "Managing models..."

# Create directory structure
create_dirs

# Train models if requested
if [ "$TRAIN" = true ]; then
    train_models
fi

# Tune models if requested
if [ "$TUNE" = true ]; then
    tune_models
fi

# Evaluate models if requested
if [ "$EVALUATE" = true ]; then
    evaluate_models
fi

# Make predictions if requested
if [ "$PREDICT" = true ]; then
    make_predictions
fi

# Optimize models if requested
if [ "$OPTIMIZE" = true ]; then
    optimize_models
fi

# Export models if requested
if [ "$EXPORT" = true ]; then
    export_models
fi

# Clean models if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean models? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_models
    fi
fi

# Backup models if requested
if [ "$BACKUP" = true ]; then
    backup_models
fi

# Restore models if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore models? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_models
    fi
fi

# Show statistics
show_stats

echo
echo "Model management completed successfully!"
