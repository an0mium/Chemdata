"""FastAPI endpoints for compound data export."""

from typing import Dict, List, Optional, Any, Set
from fastapi import APIRouter, Depends, HTTPException, Query
from fastapi.responses import JSONResponse, FileResponse, StreamingResponse
from pydantic import BaseModel, Field
import io
import csv
import json

from ...models.compound import CompoundData
from ...models.validation import ValidationError
from ..cache import Cache
from ..rate_limiter import RateLimiter
from .compounds import CompoundResponse


class ExportFormat(str):
    """Export format types."""

    JSON = "json"
    CSV = "csv"
    SDF = "sdf"
    MOL = "mol"
    SMILES = "smiles"
    INCHI = "inchi"
    EXCEL = "excel"


class ExportOptions(BaseModel):
    """Export options model."""

    format: ExportFormat = Field(..., description="Export format")
    compound_ids: Optional[List[str]] = Field(None, description="Specific compound IDs to export")
    include_fields: Optional[List[str]] = Field(None, description="Fields to include")
    exclude_fields: Optional[List[str]] = Field(None, description="Fields to exclude")
    structure_format: Optional[str] = Field(None, description="Structure output format")
    include_metadata: bool = Field(True, description="Include metadata")
    include_computed: bool = Field(True, description="Include computed properties")
    include_experimental: bool = Field(True, description="Include experimental data")
    viewport_width: Optional[int] = Field(None, description="Viewport width for responsive data")


class ExportAPI:
    """API endpoints for compound data export."""

    def __init__(self):
        """Initialize API."""
        self.router = APIRouter()
        self.cache = Cache()
        self.rate_limiter = RateLimiter()
        self._setup_routes()

    def _setup_routes(self):
        """Set up API routes."""
        self.router.add_api_route(
            "/export",
            self.export_compounds,
            methods=["POST"],
            response_class=StreamingResponse,
            tags=["export"],
        )
        self.router.add_api_route(
            "/export/formats",
            self.get_export_formats,
            methods=["GET"],
            response_model=Dict[str, Any],
            tags=["export"],
        )

    async def export_compounds(
        self,
        options: ExportOptions,
    ) -> StreamingResponse:
        """Export compounds in specified format.

        Args:
            options: Export options

        Returns:
            Streaming response with exported data

        Raises:
            HTTPException: If validation fails
        """
        # Apply rate limiting
        await self.rate_limiter.acquire()

        try:
            # Get compounds
            compounds = await self._get_compounds(options.compound_ids)

            # Format for viewport if needed
            if options.viewport_width:
                compounds = [self._format_for_viewport(c, options.viewport_width) for c in compounds]

            # Filter fields
            compounds = [
                self._filter_fields(
                    c,
                    include_fields=options.include_fields,
                    exclude_fields=options.exclude_fields,
                )
                for c in compounds
            ]

            # Format output
            if options.format == ExportFormat.JSON:
                return await self._export_json(compounds, options)
            elif options.format == ExportFormat.CSV:
                return await self._export_csv(compounds, options)
            elif options.format == ExportFormat.SDF:
                return await self._export_sdf(compounds, options)
            elif options.format == ExportFormat.MOL:
                return await self._export_mol(compounds, options)
            elif options.format == ExportFormat.SMILES:
                return await self._export_smiles(compounds, options)
            elif options.format == ExportFormat.INCHI:
                return await self._export_inchi(compounds, options)
            elif options.format == ExportFormat.EXCEL:
                return await self._export_excel(compounds, options)
            else:
                raise HTTPException(
                    status_code=400,
                    detail=f"Unsupported export format: {options.format}",
                )

        except ValidationError as e:
            raise HTTPException(status_code=400, detail=str(e))

        finally:
            self.rate_limiter.release()

    async def get_export_formats(self) -> Dict[str, Any]:
        """Get available export formats and options.

        Returns:
            Dictionary of export formats and options
        """
        return {
            "formats": [
                {
                    "id": ExportFormat.JSON,
                    "name": "JSON",
                    "description": "JSON format with full compound data",
                    "mime_type": "application/json",
                    "extension": ".json",
                },
                {
                    "id": ExportFormat.CSV,
                    "name": "CSV",
                    "description": "CSV format with tabular compound data",
                    "mime_type": "text/csv",
                    "extension": ".csv",
                },
                {
                    "id": ExportFormat.SDF,
                    "name": "SDF",
                    "description": "Structure-data file format",
                    "mime_type": "chemical/x-mdl-sdfile",
                    "extension": ".sdf",
                },
                {
                    "id": ExportFormat.MOL,
                    "name": "MOL",
                    "description": "MDL MOL file format",
                    "mime_type": "chemical/x-mdl-molfile",
                    "extension": ".mol",
                },
                {
                    "id": ExportFormat.SMILES,
                    "name": "SMILES",
                    "description": "SMILES structure format",
                    "mime_type": "text/plain",
                    "extension": ".smi",
                },
                {
                    "id": ExportFormat.INCHI,
                    "name": "InChI",
                    "description": "InChI structure format",
                    "mime_type": "text/plain",
                    "extension": ".inchi",
                },
                {
                    "id": ExportFormat.EXCEL,
                    "name": "Excel",
                    "description": "Microsoft Excel format",
                    "mime_type": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    "extension": ".xlsx",
                },
            ],
            "field_groups": [
                {
                    "id": "basic",
                    "name": "Basic Properties",
                    "fields": [
                        "id",
                        "name",
                        "smiles",
                        "inchi",
                        "molecular_weight",
                        "logp",
                        "psa",
                        "hba",
                        "hbd",
                    ],
                },
                {
                    "id": "structure",
                    "name": "Structure Data",
                    "fields": ["structure_data", "conformers", "fingerprints"],
                },
                {
                    "id": "quantum",
                    "name": "Quantum Properties",
                    "fields": ["quantum_data", "orbital_energies", "dipole"],
                },
                {
                    "id": "spectral",
                    "name": "Spectral Data",
                    "fields": ["spectral_data", "nmr_peaks", "ir_bands"],
                },
                {
                    "id": "binding",
                    "name": "Binding Data",
                    "fields": ["binding_data", "activity_data", "selectivity"],
                },
                {
                    "id": "safety",
                    "name": "Safety Data",
                    "fields": ["safety_data", "toxicity", "metabolism"],
                },
            ],
        }

    async def _get_compounds(
        self,
        compound_ids: Optional[List[str]] = None,
    ) -> List[CompoundData]:
        """Get compounds to export.

        Args:
            compound_ids: Optional list of specific compound IDs

        Returns:
            List of compound data
        """
        # TODO: Implement compound retrieval
        return []

    def _filter_fields(
        self,
        compound: CompoundData,
        include_fields: Optional[List[str]] = None,
        exclude_fields: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """Filter compound fields.

        Args:
            compound: Compound data
            include_fields: Fields to include
            exclude_fields: Fields to exclude

        Returns:
            Filtered compound data
        """
        data = compound.dict()

        if include_fields:
            return {k: v for k, v in data.items() if k in include_fields}

        if exclude_fields:
            return {k: v for k, v in data.items() if k not in exclude_fields}

        return data

    def _format_for_viewport(
        self,
        compound: CompoundData,
        viewport_width: int,
    ) -> CompoundData:
        """Format compound data for viewport.

        Args:
            compound: Compound data
            viewport_width: Viewport width

        Returns:
            Formatted compound data
        """
        # Mobile viewport
        if viewport_width < 576:
            # Reduce data for mobile
            compound.spectral_data = None
            compound.crystal_data = None
            if compound.quantum_data:
                compound.quantum_data = {k: v for k, v in compound.quantum_data.items() if k in ["energy", "dipole"]}

        # Tablet viewport
        elif viewport_width < 992:
            # Moderate data for tablet
            if compound.quantum_data:
                compound.quantum_data = {k: v for k, v in compound.quantum_data.items() if k not in ["orbital_coefficients"]}

        return compound

    async def _export_json(
        self,
        compounds: List[CompoundData],
        options: ExportOptions,
    ) -> StreamingResponse:
        """Export compounds as JSON.

        Args:
            compounds: List of compounds
            options: Export options

        Returns:
            Streaming response with JSON data
        """
        output = io.StringIO()
        json.dump([c.dict() for c in compounds], output, indent=2)
        output.seek(0)

        return StreamingResponse(
            iter([output.getvalue()]),
            media_type="application/json",
            headers={"Content-Disposition": 'attachment; filename="compounds.json"'},
        )

    async def _export_csv(
        self,
        compounds: List[CompoundData],
        options: ExportOptions,
    ) -> StreamingResponse:
        """Export compounds as CSV.

        Args:
            compounds: List of compounds
            options: Export options

        Returns:
            Streaming response with CSV data
        """
        output = io.StringIO()
        writer = csv.DictWriter(output, fieldnames=self._get_csv_fields(compounds))
        writer.writeheader()
        writer.writerows([c.dict() for c in compounds])
        output.seek(0)

        return StreamingResponse(
            iter([output.getvalue()]),
            media_type="text/csv",
            headers={"Content-Disposition": 'attachment; filename="compounds.csv"'},
        )

    async def _export_sdf(
        self,
        compounds: List[CompoundData],
        options: ExportOptions,
    ) -> StreamingResponse:
        """Export compounds as SDF.

        Args:
            compounds: List of compounds
            options: Export options

        Returns:
            Streaming response with SDF data
        """
        # TODO: Implement SDF export
        raise HTTPException(status_code=501, detail="SDF export not implemented")

    async def _export_mol(
        self,
        compounds: List[CompoundData],
        options: ExportOptions,
    ) -> StreamingResponse:
        """Export compounds as MOL files.

        Args:
            compounds: List of compounds
            options: Export options

        Returns:
            Streaming response with MOL data
        """
        # TODO: Implement MOL export
        raise HTTPException(status_code=501, detail="MOL export not implemented")

    async def _export_smiles(
        self,
        compounds: List[CompoundData],
        options: ExportOptions,
    ) -> StreamingResponse:
        """Export compounds as SMILES.

        Args:
            compounds: List of compounds
            options: Export options

        Returns:
            Streaming response with SMILES data
        """
        output = io.StringIO()
        for compound in compounds:
            if compound.structure_data and compound.structure_data.get("smiles"):
                output.write(f"{compound.structure_data['smiles']}\n")
        output.seek(0)

        return StreamingResponse(
            iter([output.getvalue()]),
            media_type="text/plain",
            headers={"Content-Disposition": 'attachment; filename="compounds.smi"'},
        )

    async def _export_inchi(
        self,
        compounds: List[CompoundData],
        options: ExportOptions,
    ) -> StreamingResponse:
        """Export compounds as InChI.

        Args:
            compounds: List of compounds
            options: Export options

        Returns:
            Streaming response with InChI data
        """
        output = io.StringIO()
        for compound in compounds:
            if compound.structure_data and compound.structure_data.get("inchi"):
                output.write(f"{compound.structure_data['inchi']}\n")
        output.seek(0)

        return StreamingResponse(
            iter([output.getvalue()]),
            media_type="text/plain",
            headers={"Content-Disposition": 'attachment; filename="compounds.inchi"'},
        )

    async def _export_excel(
        self,
        compounds: List[CompoundData],
        options: ExportOptions,
    ) -> StreamingResponse:
        """Export compounds as Excel.

        Args:
            compounds: List of compounds
            options: Export options

        Returns:
            Streaming response with Excel data
        """
        # TODO: Implement Excel export
        raise HTTPException(status_code=501, detail="Excel export not implemented")

    def _get_csv_fields(self, compounds: List[CompoundData]) -> List[str]:
        """Get CSV field names from compounds.

        Args:
            compounds: List of compounds

        Returns:
            List of field names
        """
        fields: Set[str] = set()
        for compound in compounds:
            fields.update(compound.dict().keys())
        return sorted(fields)
