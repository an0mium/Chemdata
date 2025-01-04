"""Validation logic for compound data.

This module provides the ValidationMixin class that implements validation
logic for chemical compound data.
"""

import re
from typing import List


class ValidationError(Exception):
    """Raised when compound data validation fails."""
    pass


class ValidationMixin:
    """Mixin class providing validation methods for CompoundData."""

    def validate(self) -> None:
        """Validate compound data."""
        errors = []
        errors.extend(self._validate_identifiers())
        errors.extend(self._validate_properties())
        errors.extend(self._validate_targets())
        
        if errors:
            raise ValidationError("\n".join(errors))

    def _validate_identifiers(self) -> List[str]:
        """Validate chemical identifiers."""
        errors = []
        
        # Validate CAS number
        if hasattr(self, 'cas_number') and self.cas_number and self.cas_number != "N/A":
            if not self._validate_cas_format(self.cas_number):
                errors.append(f"Invalid CAS number format: {self.cas_number}")

        # Validate SMILES
        if hasattr(self, 'smiles') and self.smiles:
            if not self._validate_smiles_format(self.smiles):
                errors.append(f"Invalid SMILES format: {self.smiles}")

        # Validate InChI
        if hasattr(self, 'inchi') and self.inchi:
            if not self._validate_inchi_format(self.inchi):
                errors.append(f"Invalid InChI format: {self.inchi}")
                
        return errors

    def _validate_properties(self) -> List[str]:
        """Validate chemical properties."""
        errors = []
        
        # Validate molecular weight
        if hasattr(self, 'molecular_weight') and self.molecular_weight < 0:
            errors.append(f"Invalid molecular weight: {self.molecular_weight}")

        # Validate LogP
        if hasattr(self, 'logp') and abs(self.logp) > 20:
            errors.append(f"Suspicious LogP value: {self.logp}")

        # Validate TPSA
        if hasattr(self, 'tpsa') and self.tpsa < 0:
            errors.append(f"Invalid TPSA value: {self.tpsa}")
            
        return errors

    def _validate_targets(self) -> List[str]:
        """Validate target data."""
        if not hasattr(self, 'targets'):
            return []
            
        errors = []
        for i, target in enumerate(self.targets):
            errors.extend(self._validate_single_target(i, target))
        return errors

    def _validate_single_target(self, index: int, target) -> List[str]:
        """Validate a single target entry."""
        errors = []
        
        # Validate affinity value
        if target.affinity_value < 0:
            errors.append(
                f"Invalid binding affinity value for target {index}: {target.affinity_value}"
            )
            
        # Validate confidence score
        if not 0 <= target.confidence <= 1:
            errors.append(
                f"Invalid confidence score for target {index}: {target.confidence}"
            )
            
        # Validate affinity type
        valid_types = {'Ki', 'IC50', 'EC50', 'Kd'}
        if target.affinity_type not in valid_types and target.affinity_type != "N/A":
            errors.append(
                f"Invalid affinity type for target {index}: {target.affinity_type}"
            )
            
        # Validate affinity unit
        valid_units = {'nM', 'uM', 'mM', 'pM'}
        if target.affinity_unit not in valid_units and target.affinity_unit != "N/A":
            errors.append(
                f"Invalid affinity unit for target {index}: {target.affinity_unit}"
            )
            
        return errors

    def _validate_cas_format(self, cas: str) -> bool:
        """
        Validate CAS number format.
        
        Args:
            cas: CAS number to validate
            
        Returns:
            True if valid, False otherwise
        """
        pattern = r'^\d{1,7}-\d{2}-\d$'
        if not re.match(pattern, cas):
            return False
            
        # Validate checksum
        numbers = cas.replace('-', '')
        check_digit = int(numbers[-1])
        numbers = numbers[:-1]
        total = sum(
            int(num) * (i + 1) 
            for i, num in enumerate(reversed(numbers))
        )
        return (total % 10) == check_digit

    def _validate_smiles_format(self, smiles: str) -> bool:
        """
        Basic validation of SMILES format.
        
        Args:
            smiles: SMILES string to validate
            
        Returns:
            True if valid format, False otherwise
        """
        # Basic format check - should contain valid element symbols
        element_pattern = r'[A-Z][a-z]?'
        if not re.search(element_pattern, smiles):
            return False
            
        # Check for balanced parentheses
        if smiles.count('(') != smiles.count(')'):
            return False
            
        # Check for balanced square brackets
        if smiles.count('[') != smiles.count(']'):
            return False
            
        return True

    def _validate_inchi_format(self, inchi: str) -> bool:
        """
        Basic validation of InChI format.
        
        Args:
            inchi: InChI string to validate
            
        Returns:
            True if valid format, False otherwise
        """
        # Should start with InChI=
        if not inchi.startswith('InChI='):
            return False
            
        # Should have at least one layer
        if not inchi.count('/') >= 1:
            return False
            
        return True
