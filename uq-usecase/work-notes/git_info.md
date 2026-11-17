# Git Workflow: Mac ↔ HPC Sync

Both Mac and HPC clones point to `origin/isatkaus/uq-usecase` on **ORNL/GridKit** (not personal fork).

---

## Mac → HPC

### 1. WIP on Mac

**Using script (recommended):**
```bash
cd /Users/isatkaus/projects/scidac/repos/gridkit
bash uq-usecase/scripts/wip_mac.sh "your commit message"
```

**Or manually:**
```bash
cd /Users/isatkaus/projects/scidac/repos/gridkit
git add uq-usecase/
git commit -m "wip: your message"
git push origin isatkaus/uq-usecase
```

### 2. Pull on HPC

```bash
cd ~/gridkit
git pull origin isatkaus/uq-usecase
```

---

## HPC → Mac

### 1. WIP on HPC

**Using script (recommended):**
```bash
cd ~/gridkit
bash uq-usecase/scripts/wip.sh "your commit message"
```

**Or manually:**
```bash
cd ~/gridkit
git add uq-usecase/
git commit -m "wip: your message"
git push origin isatkaus/uq-usecase
```

The HPC script has built-in safety: if remote has new commits with non-.md files, it warns you to pull first.

### 2. Pull on Mac

```bash
cd /Users/isatkaus/projects/scidac/repos/gridkit
git pull origin isatkaus/uq-usecase
```

If VSCode doesn't show updates immediately, reload the window (Cmd+K Cmd+R).

---

## Notes

- **Repository:** `github.com:ORNL/GridKit.git`
- **Branch:** `isatkaus/uq-usecase`
- **Main repo:** Keep separate from this branch; don't update/commit to main gridkit
- **Figures:** Restored from HPC as part of Mac workflow; always preserve when pulling
