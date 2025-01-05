# Backup Strategy

## Overview

The backup strategy needs to cover:
1. Database Backups
2. File Backups
3. Model Backups
4. Configuration Backups
5. Recovery Procedures

## Current Structure

```
backup/
└── basic_backup.sh    # Basic backup script
```

## Target Structure

```
backup/
├── database/
│   ├── full/         # Full database backups
│   └── incremental/  # Incremental backups
├── files/
│   ├── compounds/    # Compound data files
│   └── models/       # ML model files
├── config/
│   ├── system/       # System configs
│   └── app/          # App configs
└── scripts/
    ├── backup/       # Backup scripts
    └── restore/      # Restore scripts
```

## Backup Components

### 1. Database Backups

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

### 2. File Backups

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

### 3. Configuration Backups

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

### 4. Recovery Procedures

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
            
    async def restore_files(
        self,
        backup_file: Path,
        target_dir: Path
    ) -> None:
        """Restore files from backup."""
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

### Day 1: Database
1. Set up backup system
2. Configure retention
3. Test full backups
4. Test incremental

### Day 2: Files
1. Set up file backup
2. Configure archiving
3. Test compression
4. Test storage

### Day 3: Recovery
1. Set up restore
2. Test database
3. Test files
4. Test configs

### Day 4: Automation
1. Schedule backups
2. Configure monitoring
3. Set up alerts
4. Test automation

### Day 5: Documentation
1. Document procedures
2. Create runbooks
3. Test recovery
4. Train team

## Validation Steps

### 1. Backups
- [ ] Database backing up
- [ ] Files backing up
- [ ] Configs backing up
- [ ] Storage working

### 2. Recovery
- [ ] Database restores
- [ ] Files restore
- [ ] Configs restore
- [ ] Point-in-time works

### 3. Automation
- [ ] Schedules running
- [ ] Monitoring working
- [ ] Alerts firing
- [ ] Cleanup working

## Success Criteria

### 1. Reliability
- Regular backups
- Verified integrity
- Fast recovery
- No data loss

### 2. Performance
- Minimal impact
- Fast backups
- Quick recovery
- Good compression

### 3. Usability
- Easy recovery
- Clear procedures
- Good monitoring
- Fast verification

## Next Steps

1. Set up infrastructure
2. Configure backups
3. Test recovery
4. Document procedures
5. Train team
6. Monitor system
