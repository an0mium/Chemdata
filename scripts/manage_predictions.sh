#!/bin/bash
# Script to manage ML predictions, inference, and analysis

# Exit on error
set -e

# Default values
DATA_DIR="data/predictions"
MODEL_DIR="models"
PRED_DIR="predictions"
CACHE_DIR="cache/predictions"
LOG_DIR="logs"
PREDICTION_TYPES="activity,toxicity,abuse,psychoactive,nootropic,bbb,binding,affinity"
MODEL_STAGES="base,tuned,ensemble"
MODEL_MODES="development,production,testing"
PRED_MODES="single,batch,realtime"
PRED_STAGES="raw,processed,validated"
DATASET_TYPES="bindingdb,chembl,pubchem,swiss,community,literature"
BATCH_SIZE=32
THRESHOLD=0.5
CONFIDENCE=0.7
TRAIN=false
PREDICT=false
ANALYZE=false
VALIDATE=false
EVALUATE=false
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
        --pred-dir)
            PRED_DIR="$2"
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
        --prediction-types)
            PREDICTION_TYPES="$2"
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
        --pred-modes)
            PRED_MODES="$2"
            shift 2
            ;;
        --pred-stages)
            PRED_STAGES="$2"
            shift 2
            ;;
        --dataset-types)
            DATASET_TYPES="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --threshold)
            THRESHOLD="$2"
            shift 2
            ;;
        --confidence)
            CONFIDENCE="$2"
            shift 2
            ;;
        --train)
            TRAIN=true
            shift
            ;;
        --predict)
            PREDICT=true
            shift
            ;;
        --analyze)
            ANALYZE=true
            shift
            ;;
        --validate)
            VALIDATE=true
            shift
            ;;
        --evaluate)
            EVALUATE=true
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
    for type in ${PREDICTION_TYPES//,/ }; do
        for stage in ${MODEL_STAGES//,/ }; do
            for mode in ${MODEL_MODES//,/ }; do
                mkdir -p "$MODEL_DIR/$type/$stage/$mode"/{checkpoints,configs,metrics,reports}
            done
        done
    done
    
    # Prediction directories
    for type in ${PREDICTION_TYPES//,/ }; do
        for mode in ${PRED_MODES//,/ }; do
            for stage in ${PRED_STAGES//,/ }; do
                mkdir -p "$PRED_DIR/$type/$mode/$stage"/{data,reports}
            done
        done
    done
    
    # Data directories
    for type in ${DATASET_TYPES//,/ }; do
        mkdir -p "$DATA_DIR/$type"/{raw,processed,analyzed,validated,reports}
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
        
        for type in ${PREDICTION_TYPES//,/ }; do
            echo "Training $type model..."
            
            for stage in ${MODEL_STAGES//,/ }; do
                echo "Training $stage model..."
                
                # Build command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli train-model"
                CMD="$CMD --data-dir $DATA_DIR"
                CMD="$CMD --model-dir $MODEL_DIR/$type/$stage/$mode"
                CMD="$CMD --cache-dir $CACHE_DIR"
                CMD="$CMD --type $type"
                CMD="$CMD --stage $stage"
                CMD="$CMD --mode $mode"
                CMD="$CMD --batch-size $BATCH_SIZE"
                CMD="$CMD --threshold $THRESHOLD"
                CMD="$CMD --confidence $CONFIDENCE"
                CMD="$CMD --log-dir $LOG_DIR"
                
                if [ "$GPU" = true ]; then
                    CMD="$CMD --gpu"
                fi
                
                # Run command
                echo "Running: $CMD"
                $CMD || true  # Continue even if training fails
                
                # Generate training report
                echo "Generating training report..."
                python -m binding_data_processor.processors.psychopharm.predictors.cli generate-training-report \
                    --model-dir "$MODEL_DIR/$type/$stage/$mode" \
                    --output "$MODEL_DIR/$type/$stage/$mode/reports/training.html"
            done
        done
    done
}

# Function to make predictions
make_predictions() {
    echo "Making predictions..."
    
    for pred_mode in ${PRED_MODES//,/ }; do
        echo "Making predictions in $pred_mode mode..."
        
        for type in ${PREDICTION_TYPES//,/ }; do
            echo "Making $type predictions..."
            
            for stage in ${MODEL_STAGES//,/ }; do
                echo "Making predictions with $stage model..."
                
                for model_mode in ${MODEL_MODES//,/ }; do
                    echo "Using model from $model_mode mode..."
                    
                    # Predict on each dataset type
                    for dataset in ${DATASET_TYPES//,/ }; do
                        echo "Predicting on $dataset dataset..."
                        
                        # Build prediction command
                        CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli predict"
                        CMD="$CMD --data-dir $DATA_DIR/$dataset"
                        CMD="$CMD --model-dir $MODEL_DIR/$type/$stage/$model_mode"
                        CMD="$CMD --pred-dir $PRED_DIR/$type/$pred_mode/raw"
                        CMD="$CMD --cache-dir $CACHE_DIR"
                        CMD="$CMD --type $type"
                        CMD="$CMD --stage $stage"
                        CMD="$CMD --model-mode $model_mode"
                        CMD="$CMD --pred-mode $pred_mode"
                        CMD="$CMD --dataset $dataset"
                        CMD="$CMD --batch-size $BATCH_SIZE"
                        CMD="$CMD --threshold $THRESHOLD"
                        CMD="$CMD --confidence $CONFIDENCE"
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
                            --pred-dir "$PRED_DIR/$type/$pred_mode/raw" \
                            --dataset "$dataset" \
                            --output "$PRED_DIR/$type/$pred_mode/raw/reports/prediction_${dataset}.html"
                    done
                done
            done
        done
    done
}

# Function to analyze predictions
analyze_predictions() {
    echo "Analyzing predictions..."
    
    for pred_mode in ${PRED_MODES//,/ }; do
        echo "Analyzing predictions in $pred_mode mode..."
        
        for type in ${PREDICTION_TYPES//,/ }; do
            echo "Analyzing $type predictions..."
            
            # Build analysis command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli analyze-predictions"
            CMD="$CMD --pred-dir $PRED_DIR/$type/$pred_mode"
            CMD="$CMD --type $type"
            CMD="$CMD --mode $pred_mode"
            CMD="$CMD --threshold $THRESHOLD"
            CMD="$CMD --confidence $CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run analysis
            echo "Running: $CMD"
            $CMD || true  # Continue even if analysis fails
            
            # Generate analysis report
            echo "Generating analysis report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-analysis-report \
                --pred-dir "$PRED_DIR/$type/$pred_mode" \
                --output "$PRED_DIR/$type/$pred_mode/processed/reports/analysis.html"
        done
    done
}

# Function to validate predictions
validate_predictions() {
    echo "Validating predictions..."
    
    for pred_mode in ${PRED_MODES//,/ }; do
        echo "Validating predictions in $pred_mode mode..."
        
        for type in ${PREDICTION_TYPES//,/ }; do
            echo "Validating $type predictions..."
            
            # Build validation command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli validate-predictions"
            CMD="$CMD --pred-dir $PRED_DIR/$type/$pred_mode"
            CMD="$CMD --type $type"
            CMD="$CMD --mode $pred_mode"
            CMD="$CMD --threshold $THRESHOLD"
            CMD="$CMD --confidence $CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run validation
            echo "Running: $CMD"
            $CMD || true  # Continue even if validation fails
            
            # Generate validation report
            echo "Generating validation report..."
            python -m binding_data_processor.processors.psychopharm.predictors.cli generate-validation-report \
                --pred-dir "$PRED_DIR/$type/$pred_mode" \
                --output "$PRED_DIR/$type/$pred_mode/validated/reports/validation.html"
        done
    done
}

# Function to evaluate models
evaluate_models() {
    echo "Evaluating models..."
    
    for mode in ${MODEL_MODES//,/ }; do
        echo "Evaluating models in $mode mode..."
        
        for type in ${PREDICTION_TYPES//,/ }; do
            echo "Evaluating $type model..."
            
            for stage in ${MODEL_STAGES//,/ }; do
                echo "Evaluating $stage model..."
                
                # Build evaluation command
                CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli evaluate-model"
                CMD="$CMD --data-dir $DATA_DIR"
                CMD="$CMD --model-dir $MODEL_DIR/$type/$stage/$mode"
                CMD="$CMD --type $type"
                CMD="$CMD --stage $stage"
                CMD="$CMD --mode $mode"
                CMD="$CMD --batch-size $BATCH_SIZE"
                CMD="$CMD --threshold $THRESHOLD"
                CMD="$CMD --confidence $CONFIDENCE"
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
                    --output "$MODEL_DIR/$type/$stage/$mode/reports/evaluation.html"
            done
        done
    done
}

# Function to export predictions
export_predictions() {
    echo "Exporting predictions..."
    
    for pred_mode in ${PRED_MODES//,/ }; do
        echo "Exporting predictions in $pred_mode mode..."
        
        for type in ${PREDICTION_TYPES//,/ }; do
            echo "Exporting $type predictions..."
            
            # Build export command
            CMD="python -m binding_data_processor.processors.psychopharm.predictors.cli export-predictions"
            CMD="$CMD --pred-dir $PRED_DIR/$type/$pred_mode"
            CMD="$CMD --type $type"
            CMD="$CMD --mode $pred_mode"
            CMD="$CMD --threshold $THRESHOLD"
            CMD="$CMD --confidence $CONFIDENCE"
            CMD="$CMD --log-dir $LOG_DIR"
            
            # Run export
            echo "Running: $CMD"
            $CMD || true  # Continue even if export fails
        done
    done
}

# Function to clean predictions
clean_predictions() {
    echo "Cleaning predictions..."
    
    # Clean model directories
    rm -rf "$MODEL_DIR"/*
    
    # Clean prediction directories
    rm -rf "$PRED_DIR"/*
    
    # Clean data directories
    rm -rf "$DATA_DIR"/*
    
    # Clean cache
    rm -rf "$CACHE_DIR"/*
}

# Function to backup predictions
backup_predictions() {
    echo "Backing up predictions..."
    
    # Create backup directory
    BACKUP_DIR="backups/predictions_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$BACKUP_DIR"
    
    # Copy directories
    cp -r "$MODEL_DIR" "$BACKUP_DIR/"
    cp -r "$PRED_DIR" "$BACKUP_DIR/"
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

# Function to restore predictions
restore_predictions() {
    echo "Restoring predictions..."
    
    # Find latest backup
    local LATEST_BACKUP=""
    if command_exists gpg; then
        LATEST_BACKUP=$(ls -t backups/predictions_*.tar.gz.gpg 2>/dev/null | head -n1)
        if [ -n "$LATEST_BACKUP" ]; then
            echo "Decrypting backup: $LATEST_BACKUP"
            gpg --decrypt "$LATEST_BACKUP" | tar -xz
            BACKUP_DIR="${LATEST_BACKUP%.tar.gz.gpg}"
        fi
    fi
    
    if [ -z "$LATEST_BACKUP" ]; then
        LATEST_BACKUP=$(ls -t backups/predictions_*.tar.gz 2>/dev/null | head -n1)
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
        rm -rf "$MODEL_DIR" "$PRED_DIR" "$DATA_DIR"
    fi
    
    cp -r "$BACKUP_DIR/models" ./
    cp -r "$BACKUP_DIR/predictions" ./
    cp -r "$BACKUP_DIR/data/predictions" ./data/
    
    # Clean up
    rm -rf "$BACKUP_DIR"
    
    echo "Predictions restored from: $BACKUP_DIR"
}

# Function to show prediction statistics
show_stats() {
    echo "Prediction statistics:"
    echo
    
    echo "Model files:"
    for mode in ${MODEL_MODES//,/ }; do
        echo "$mode mode:"
        for type in ${PREDICTION_TYPES//,/ }; do
            echo "  $type model:"
            
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
                
                # Show reports
                if [ -f "$MODEL_DIR/$type/$stage/$mode/reports/training.html" ]; then
                    echo "      Training report: $MODEL_DIR/$type/$stage/$mode/reports/training.html"
                fi
                if [ -f "$MODEL_DIR/$type/$stage/$mode/reports/evaluation.html" ]; then
                    echo "      Evaluation report: $MODEL_DIR/$type/$stage/$mode/reports/evaluation.html"
                fi
                echo
            done
        done
    done
    
    echo "Prediction files:"
    for pred_mode in ${PRED_MODES//,/ }; do
        echo "$pred_mode mode:"
        for type in ${PREDICTION_TYPES//,/ }; do
            echo "  $type predictions:"
            
            for stage in ${PRED_STAGES//,/ }; do
                echo "    $stage stage:"
                echo "      Files: $(find "$PRED_DIR/$type/$pred_mode/$stage/data" -type f | wc -l) files"
                echo "      Size: $(du -sh "$PRED_DIR/$type/$pred_mode/$stage/data" | cut -f1)"
                
                # Count predictions by type
                if [ -f "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json" ]; then
                    case "$type" in
                        activity)
                            echo "      Activities: $(jq '.predictions | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            echo "      Targets: $(jq '.targets | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            if [ "$stage" = "processed" ]; then
                                echo "      High confidence: $(jq --arg t "$CONFIDENCE" '.predictions[] | select(.confidence >= ($t | tonumber)) | .id' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json" | wc -l)"
                            fi
                            ;;
                        toxicity)
                            echo "      Toxicity: $(jq '.predictions | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            echo "      Effects: $(jq '.effects | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            if [ "$stage" = "processed" ]; then
                                echo "      High risk: $(jq --arg t "$CONFIDENCE" '.predictions[] | select(.risk >= ($t | tonumber)) | .id' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json" | wc -l)"
                            fi
                            ;;
                        abuse)
                            echo "      Abuse: $(jq '.predictions | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            echo "      Mechanisms: $(jq '.mechanisms | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            if [ "$stage" = "processed" ]; then
                                echo "      High potential: $(jq --arg t "$CONFIDENCE" '.predictions[] | select(.potential >= ($t | tonumber)) | .id' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json" | wc -l)"
                            fi
                            ;;
                        psychoactive)
                            echo "      Effects: $(jq '.predictions | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            echo "      Classes: $(jq '.classes | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            if [ "$stage" = "processed" ]; then
                                echo "      High activity: $(jq --arg t "$CONFIDENCE" '.predictions[] | select(.activity >= ($t | tonumber)) | .id' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json" | wc -l)"
                            fi
                            ;;
                        nootropic)
                            echo "      Effects: $(jq '.predictions | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            echo "      Mechanisms: $(jq '.mechanisms | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            if [ "$stage" = "processed" ]; then
                                echo "      High efficacy: $(jq --arg t "$CONFIDENCE" '.predictions[] | select(.efficacy >= ($t | tonumber)) | .id' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json" | wc -l)"
                            fi
                            ;;
                        bbb)
                            echo "      BBB: $(jq '.predictions | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            echo "      Properties: $(jq '.properties | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            if [ "$stage" = "processed" ]; then
                                echo "      High penetration: $(jq --arg t "$CONFIDENCE" '.predictions[] | select(.penetration >= ($t | tonumber)) | .id' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json" | wc -l)"
                            fi
                            ;;
                        binding)
                            echo "      Binding: $(jq '.predictions | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            echo "      Targets: $(jq '.targets | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            if [ "$stage" = "processed" ]; then
                                echo "      High affinity: $(jq --arg t "$CONFIDENCE" '.predictions[] | select(.affinity >= ($t | tonumber)) | .id' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json" | wc -l)"
                            fi
                            ;;
                        affinity)
                            echo "      Affinity: $(jq '.predictions | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            echo "      Targets: $(jq '.targets | length' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json")"
                            if [ "$stage" = "processed" ]; then
                                echo "      High affinity: $(jq --arg t "$CONFIDENCE" '.predictions[] | select(.affinity >= ($t | tonumber)) | .id' "$PRED_DIR/$type/$pred_mode/$stage/data/predictions.json" | wc -l)"
                            fi
                            ;;
                    esac
                    
                    if [ "$stage" = "validated" ] && [ -f "$PRED_DIR/$type/$pred_mode/$stage/data/validation.json" ]; then
                        echo "      Validation score: $(jq '.validation_score' "$PRED_DIR/$type/$pred_mode/$stage/data/validation.json")"
                    fi
                fi
                
                # Show reports
                if [ -f "$PRED_DIR/$type/$pred_mode/$stage/reports/predictions.html" ]; then
                    echo "      Prediction report: $PRED_DIR/$type/$pred_mode/$stage/reports/predictions.html"
                fi
                if [ -f "$PRED_DIR/$type/$pred_mode/$stage/reports/validation.html" ]; then
                    echo "      Validation report: $PRED_DIR/$type/$pred_mode/$stage/reports/validation.html"
                fi
            done
            echo
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
        
        # Analyzed data
        echo "  Analyzed data:"
        echo "    Files: $(find "$DATA_DIR/$type/analyzed" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DATA_DIR/$type/analyzed" | cut -f1)"
        
        # Validated data
        echo "  Validated data:"
        echo "    Files: $(find "$DATA_DIR/$type/validated" -type f | wc -l) files"
        echo "    Size: $(du -sh "$DATA_DIR/$type/validated" | cut -f1)"
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
echo "Managing predictions..."

# Create directory structure
create_dirs

# Train models if requested
if [ "$TRAIN" = true ]; then
    train_models
fi

# Make predictions if requested
if [ "$PREDICT" = true ]; then
    make_predictions
fi

# Analyze predictions if requested
if [ "$ANALYZE" = true ]; then
    analyze_predictions
fi

# Validate predictions if requested
if [ "$VALIDATE" = true ]; then
    validate_predictions
fi

# Evaluate models if requested
if [ "$EVALUATE" = true ]; then
    evaluate_models
fi

# Export predictions if requested
if [ "$EXPORT" = true ]; then
    export_predictions
fi

# Clean predictions if requested
if [ "$CLEAN" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to clean predictions? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        clean_predictions
    fi
fi

# Backup predictions if requested
if [ "$BACKUP" = true ]; then
    backup_predictions
fi

# Restore predictions if requested
if [ "$RESTORE" = true ]; then
    if [ "$FORCE" = true ] || read -p "Are you sure you want to restore predictions? [y/N] " -n 1 -r && [[ $REPLY =~ ^[Yy]$ ]]; then
        echo
        restore_predictions
    fi
fi

# Show statistics
show_stats

echo
echo "Prediction management completed successfully!"
