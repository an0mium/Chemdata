# Backup Strategy

## Overview

The backup strategy needs to cover:
1. Database Backups (Priority)
2. Web Assets Backups (Priority)
3. File Backups
4. Model Backups
5. Configuration Backups
6. Recovery Procedures

## Current Structure

```
backup/
└── basic_backup.sh    # Basic backup script
```

## Target Structure

```
backup/
├── database/         # Database backups (Priority)
│   ├── full/        # Full database backups
│   └── incremental/ # Incremental backups
├── web/             # Web assets (Priority)
│   ├── static/      # Static assets (CSS, JS, images)
│   ├── templates/   # HTML templates
│   └── responsive/  # Responsive design assets
├── files/
│   ├── compounds/   # Compound data files
│   └── models/      # ML model files
├── config/
│   ├── system/      # System configs
│   └── app/         # App configs
└── scripts/
    ├── backup/      # Backup scripts
    └── restore/     # Restore scripts
```

## Backup Components

### 1. Database Backups (Priority)

```python
# In backup/database/manager.py
class DatabaseBackup:
    """Database backup management."""
    def __init__(self):
        self.config = BackupConfig()
        self.storage = StorageManager()
        
    async def create_full_backup(self) -> Path:
        """Create full database backup."""
        # Get connection
        async with self.get_db_connection() as conn:
            # Create backup
            backup_file = self.config.backup_dir / f"db_full_{timestamp()}.sql"
            
            # Dump database
            process = await asyncio.create_subprocess_exec(
                "pg_dump",
                "-h", self.config.db_host,
                "-U", self.config.db_user,
                "-F", "c",  # Custom format
                "-f", str(backup_file),
                self.config.db_name,
                env={"PGPASSWORD": self.config.db_password}
            )
            
            # Wait for completion
            await process.wait()
            
            # Upload to storage
            await self.storage.upload_file(
                backup_file,
                f"database/full/{backup_file.name}"
            )
            
            return backup_file
            
    async def create_incremental_backup(self) -> Path:
        """Create incremental backup."""
        # Get last full backup
        last_full = await self.get_last_full_backup()
        
        # Get changes since last backup
        async with self.get_db_connection() as conn:
            backup_file = self.config.backup_dir / f"db_inc_{timestamp()}.sql"
            
            # Dump changes
            process = await asyncio.create_subprocess_exec(
                "pg_dump",
                "-h", self.config.db_host,
                "-U", self.config.db_user,
                "--since", last_full.timestamp,
                "-F", "c",
                "-f", str(backup_file),
                self.config.db_name,
                env={"PGPASSWORD": self.config.db_password}
            )
            
            await process.wait()
            
            # Upload to storage
            await self.storage.upload_file(
                backup_file,
                f"database/incremental/{backup_file.name}"
            )
            
            return backup_file
```

### 2. Web Assets Backups (Priority)

```python
# In backup/web/manager.py
class WebAssetsBackup:
    """Web assets backup management."""
    def __init__(self):
        self.config = BackupConfig()
        self.storage = StorageManager()
        
    async def backup_static_assets(self) -> Path:
        """Backup static web assets."""
        # Get static files
        static_dir = self.config.web_dir / "static"
        
        # Create archive
        archive = self.config.backup_dir / f"web_static_{timestamp()}.tar.gz"
        
        # Compress files
        process = await asyncio.create_subprocess_exec(
            "tar",
            "-czf", str(archive),
            "-C", str(static_dir),
            "."
        )
        
        await process.wait()
        
        # Upload to storage
        await self.storage.upload_file(
            archive,
            f"web/static/{archive.name}"
        )
        
        return archive
        
    async def backup_templates(self) -> Path:
        """Backup HTML templates."""
        # Get template files
        template_dir = self.config.web_dir / "templates"
        
        # Create archive
        archive = self.config.backup_dir / f"web_templates_{timestamp()}.tar.gz"
        
        # Compress files
        process = await asyncio.create_subprocess_exec(
            "tar",
            "-czf", str(archive),
            "-C", str(template_dir),
            "."
        )
        
        await process.wait()
        
        # Upload to storage
        await self.storage.upload_file(
            archive,
            f"web/templates/{archive.name}"
        )
        
        return archive
        
    async def backup_responsive_assets(self) -> Path:
        """Backup responsive design assets."""
        # Get responsive design files
        responsive_dir = self.config.web_dir / "static" / "css" / "responsive"
        
        # Create archive
        archive = self.config.backup_dir / f"web_responsive_{timestamp()}.tar.gz"
        
        # Compress files
        process = await asyncio.create_subprocess_exec(
            "tar",
            "-czf", str(archive),
            "-C", str(responsive_dir),
            "."
        )
        
        await process.wait()
        
        # Upload to storage
        await self.storage.upload_file(
            archive,
            f"web/responsive/{archive.name}"
        )
        
        return archive
```

### 3. File & Model Backups

```python
# In backup/files/manager.py
class FileBackup:
    """File backup management."""
    def __init__(self):
        self.config = BackupConfig()
        self.storage = StorageManager()
        
    async def backup_compounds(self) -> List[Path]:
        """Backup compound data files."""
        # Get compound files
        compound_dir = self.config.data_dir / "compounds"
        files = list(compound_dir.glob("**/*.tsv"))
        
        # Create archive
        archive = self.config.backup_dir / f"compounds_{timestamp()}.tar.gz"
        
        # Compress files
        process = await asyncio.create_subprocess_exec(
            "tar",
            "-czf", str(archive),
            *[str(f) for f in files],
            cwd=str(compound_dir)
        )
        
        await process.wait()
        
        # Upload to storage
        await self.storage.upload_file(
            archive,
            f"files/compounds/{archive.name}"
        )
        
        return archive
        
    async def backup_models(self) -> List[Path]:
        """Backup ML model files."""
        # Get model files
        model_dir = self.config.model_dir
        files = list(model_dir.glob("**/*.pkl"))
        
        # Create archive
        archive = self.config.backup_dir / f"models_{timestamp()}.tar.gz"
        
        # Compress files
        process = await asyncio.create_subprocess_exec(
            "tar",
            "-czf", str(archive),
            *[str(f) for f in files],
            cwd=str(model_dir)
        )
        
        await process.wait()
        
        # Upload to storage
        await self.storage.upload_file(
            archive,
            f"files/models/{archive.name}"
        )
        
        return archive
```

### 4. Configuration Backups

```python
# In backup/config/manager.py
class ConfigBackup:
    """Configuration backup management."""
    def __init__(self):
        self.config = BackupConfig()
        self.storage = StorageManager()
        
    async def backup_configs(self) -> Path:
        """Backup all configuration files."""
        # Get config files
        config_files = [
            *self.config.config_dir.glob("*.yaml"),
            *self.config.config_dir.glob("*.env"),
            *self.config.config_dir.glob("*.json")
        ]
        
        # Create archive
        archive = self.config.backup_dir / f"config_{timestamp()}.tar.gz"
        
        # Compress files
        process = await asyncio.create_subprocess_exec(
            "tar",
            "-czf", str(archive),
            *[str(f) for f in config_files],
            cwd=str(self.config.config_dir)
        )
        
        await process.wait()
        
        # Upload to storage
        await self.storage.upload_file(
            archive,
            f"config/system/{archive.name}"
        )
        
        return archive
```

### 5. Recovery Procedures

```python
# In backup/restore/manager.py
class RestoreManager:
    """Backup restoration management."""
    def __init__(self):
        self.config = BackupConfig()
        self.storage = StorageManager()
        
    async def restore_database(
        self,
        backup_file: Path,
        target_time: Optional[datetime] = None
    ) -> None:
        """Restore database from backup."""
        # Download backup
        local_file = await self.storage.download_file(backup_file)
        
        # Restore database
        process = await asyncio.create_subprocess_exec(
            "pg_restore",
            "-h", self.config.db_host,
            "-U", self.config.db_user,
            "-d", self.config.db_name,
            "-F", "c",
            str(local_file),
            env={"PGPASSWORD": self.config.db_password}
        )
        
        await process.wait()
        
        # Apply incremental if needed
        if target_time:
            await self.apply_incremental_backups(target_time)
            
    async def restore_web_assets(
        self,
        backup_file: Path,
        target_dir: Path
    ) -> None:
        """Restore web assets from backup."""
        # Download backup
        local_file = await self.storage.download_file(backup_file)
        
        # Extract files
        process = await asyncio.create_subprocess_exec(
            "tar",
            "-xzf", str(local_file),
            "-C", str(target_dir)
        )
        
        await process.wait()
```

## Implementation Steps

### Day 1: Priority Backups
1. Set up database backup system
2. Configure web assets backup
3. Test database backups
4. Test web assets backups
5. Document procedures

### Day 2: File & Model Backups
1. Set up file backup
2. Configure model backup
3. Test compression
4. Test storage
5. Document procedures

### Day 3: Configuration & Recovery
1. Set up config backup
2. Configure recovery
3. Test restore
4. Document procedures

### Day 4: Integration & Automation
1. Schedule backups
2. Configure monitoring
3. Set up alerts
4. Test automation

### Day 5: Testing & Documentation
1. Test all backups
2. Test all recovery
3. Update docs
4. Train team

## Validation Steps

### 1. Database Backups (Priority)
- [ ] Full backups working
- [ ] Incremental working
- [ ] Recovery tested
- [ ] Point-in-time works

### 2. Web Assets (Priority)
- [ ] Static assets backing up
- [ ] Templates backing up
- [ ] Responsive assets backing up
- [ ] Recovery tested

### 3. File & Model Backups
- [ ] Files backing up
- [ ] Models backing up
- [ ] Recovery tested
- [ ] Storage working

### 4. Configuration & Recovery
- [ ] Configs backing up
- [ ] Recovery working
- [ ] Procedures documented
- [ ] Team trained

### 5. Integration & Automation
- [ ] Schedules running
- [ ] Monitoring working
- [ ] Alerts firing
- [ ] Cleanup working

## Success Criteria

### 1. Database Health (Priority)
- Regular backups
- Point-in-time recovery
- Data integrity
- Fast recovery
- Good compression

### 2. Web Assets (Priority)
- Regular backups
- Fast recovery
- Asset integrity
- Style preservation
- Responsive layouts

### 3. File Management
- Regular backups
- Fast recovery
- Data integrity
- Good compression
- Easy access

### 4. Recovery Speed
- Quick restore
- Data integrity
- Asset preservation
- Easy rollback
- Clear procedures

### 5. Integration Quality
- Automated backups
- Good monitoring
- Clear alerts
- Easy management
- Fast verification

### 6. Reliability
- Regular backups
- Verified integrity
- Fast recovery
- No data loss
- Clear procedures

### 7. Performance
- Minimal impact
- Fast backups
- Quick recovery
- Good compression
- Efficient storage

### 8. Usability
- Easy recovery
- Clear procedures
- Good monitoring
- Fast verification
- Simple management

## Next Steps

1. Set up database backups
2. Configure web assets backups
3. Add file backups
4. Set up config backups
5. Test recovery
6. Configure automation
7. Set up monitoring
8. Document procedures
9. Train team
