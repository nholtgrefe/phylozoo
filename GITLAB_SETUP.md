# GitLab Setup Guide

This guide will help you set up Git and connect your repository to GitLab.

## Prerequisites

1. **Git installed** - Check with: `git --version`
2. **GitLab account** - You need access to a GitLab instance (e.g., gitlab.tudelft.nl)

## Initial Setup

### Step 0: Create Repository on GitLab (Do This First!)

**Yes, create the repository on GitLab first!** Here's how:

1. **Log in to GitLab** (e.g., https://gitlab.tudelft.nl)
2. **Click "New project"** or the "+" button
3. **Choose "Create blank project"**
4. **Fill in the details:**
   - Project name: `phylozoo` (or your preferred name)
   - Project slug: (auto-generated, usually fine)
   - Visibility: Choose Private or Internal (Public if you want it open)
   - **IMPORTANT: Do NOT initialize with a README, .gitignore, or license** (since you already have these files)
5. **Click "Create project"**

6. **Copy the repository URL** - GitLab will show you the URL. It will look like:
   - HTTPS: `https://gitlab.tudelft.nl/YOUR_USERNAME/phylozoo.git`
   - SSH: `git@gitlab.tudelft.nl:YOUR_USERNAME/phylozoo.git`

You'll need this URL in the next steps!

### Step 1: Initialize Git Repository

If your project isn't already a Git repository:

```bash
cd /home/nholtgreve/Documents/Code/phylozoo
git init
```

### Step 2: Configure Git (if not already done)

Set your name and email (use your GitLab email):

```bash
git config --global user.name "N. Holtgrefe"
git config --global user.email "n.a.l.holtgrefe@tudelft.nl"
```

Or set it only for this repository:

```bash
git config user.name "N. Holtgrefe"
git config user.email "n.a.l.holtgrefe@tudelft.nl"
```

### Step 3: Add GitLab Remote

Connect your local repository to GitLab. Replace `YOUR_USERNAME` and `YOUR_REPO_NAME` with your actual values:

**For HTTPS (recommended for beginners):**
```bash
git remote add origin https://gitlab.tudelft.nl/YOUR_USERNAME/YOUR_REPO_NAME.git
```

**For SSH (if you have SSH keys set up):**
```bash
git remote add origin git@gitlab.tudelft.nl:YOUR_USERNAME/YOUR_REPO_NAME.git
```

To check your remote:
```bash
git remote -v
```

### Step 4: Stage and Commit Your Files

Add all files to Git:

```bash
git add .
```

Or add specific files:
```bash
git add README.md pyproject.toml src/ tests/
```

Create your first commit:

```bash
git commit -m "Initial commit: PhyloZoo package setup"
```

### Step 5: Push to GitLab

Push your code to GitLab:

```bash
git push -u origin main
```

If your default branch is `master` instead of `main`:
```bash
git branch -M main  # Rename branch to main
git push -u origin main
```

## Daily Workflow

### Making Changes

1. **Check status** - See what files have changed:
   ```bash
   git status
   ```

2. **Stage changes** - Add files you want to commit:
   ```bash
   git add <file1> <file2>
   # Or add all changes:
   git add .
   ```

3. **Commit changes** - Save your changes with a descriptive message:
   ```bash
   git commit -m "Description of what you changed"
   ```

4. **Push to GitLab** - Upload your commits:
   ```bash
   git push
   ```

### Pulling Changes

If working with others or on multiple machines, pull latest changes:

```bash
git pull
```

### Viewing History

```bash
git log                    # View commit history
git log --oneline          # Compact view
git log --graph --oneline  # With branch visualization
```

## Branching (Recommended for Collaboration)

### Create a New Branch

```bash
git checkout -b feature/my-new-feature
# Or with newer Git:
git switch -c feature/my-new-feature
```

### Switch Between Branches

```bash
git checkout main
# Or:
git switch main
```

### Push a Branch to GitLab

```bash
git push -u origin feature/my-new-feature
```

### Merge Branches

```bash
git checkout main
git merge feature/my-new-feature
```

## Common GitLab Workflows

### Creating Merge Requests (Pull Requests)

1. Push your branch to GitLab:
   ```bash
   git push -u origin feature/my-feature
   ```

2. Go to GitLab web interface
3. Click "Create merge request" when prompted
4. Fill in the merge request details
5. Request review from team members

### Cloning a Repository

If you need to clone this repository on another machine:

```bash
git clone https://gitlab.tudelft.nl/YOUR_USERNAME/YOUR_REPO_NAME.git
cd YOUR_REPO_NAME
```

### Updating from Remote

```bash
git fetch origin          # Download changes without merging
git pull origin main      # Download and merge changes
```

## Useful Git Commands

```bash
# See what changed
git diff

# See staged changes
git diff --staged

# Undo changes to a file (before staging)
git restore <file>

# Unstage a file
git restore --staged <file>

# View remote repositories
git remote -v

# Update remote URL
git remote set-url origin <new-url>

# Create a tag (for releases)
git tag -a v1.0.0 -m "Version 1.0.0"
git push origin v1.0.0
```

## Authentication

### HTTPS Authentication

GitLab will prompt for credentials. For better security, use a Personal Access Token:

1. Go to GitLab → Settings → Access Tokens
2. Create a token with `read_repository` and `write_repository` scopes
3. Use the token as your password when prompted

### SSH Authentication (Recommended)

1. Generate SSH key (if you don't have one):
   ```bash
   ssh-keygen -t ed25519 -C "n.a.l.holtgrefe@tudelft.nl"
   ```

2. Add to SSH agent:
   ```bash
   eval "$(ssh-agent -s)"
   ssh-add ~/.ssh/id_ed25519
   ```

3. Copy public key:
   ```bash
   cat ~/.ssh/id_ed25519.pub
   ```

4. Add to GitLab: Settings → SSH Keys → Add key

## Troubleshooting

### "Repository not found" error
- Check the repository URL
- Verify you have access to the repository
- Check authentication

### "Permission denied" error
- Check SSH keys are set up correctly
- Verify your GitLab account has access
- Try using HTTPS instead

### Merge conflicts
```bash
# See conflicted files
git status

# Resolve conflicts in files, then:
git add <resolved-files>
git commit -m "Resolve merge conflicts"
```

## Best Practices

1. **Commit often** - Small, frequent commits are better than large ones
2. **Write clear commit messages** - Describe what and why, not just what
3. **Use branches** - Keep main branch stable, work on feature branches
4. **Pull before push** - Always pull latest changes before pushing
5. **Review before merging** - Use merge requests for code review

## Quick Reference

```bash
# Daily workflow
git status                    # Check what changed
git add .                     # Stage all changes
git commit -m "Message"       # Commit changes
git push                      # Push to GitLab

# Working with branches
git checkout -b new-branch    # Create and switch to branch
git checkout main             # Switch to main branch
git merge branch-name         # Merge branch into current

# Getting updates
git pull                      # Pull latest changes
git fetch                     # Download without merging
```

