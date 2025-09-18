# Security in BioVault

## SyftBox Permission System

BioVault leverages SyftBox's powerful Access Control List (ACL) system to protect sensitive biological data. This document explains how permissions work and how BioVault ensures data security.

## How SyftBox Protects Data

SyftBox uses a hierarchical, file-based permission system with YAML configuration files (`syft.pub.yaml`) placed in directories throughout the datasite structure. Each file controls access to its directory and subdirectories.

### Key Security Features

1. **Owner-Based Access Control**: Datasite owners have implicit full access to their files
2. **Hierarchical Rules**: Permissions cascade from parent to child directories
3. **Terminal Nodes**: Can prevent subdirectories from overriding parent permissions
4. **Template-Based Patterns**: Dynamic access control based on user identity
5. **Granular Permissions**: Separate controls for read, write, and admin operations

## Permission Test Scenarios

BioVault includes comprehensive permission tests to ensure data security. Here's how the system protects data through various scenarios:

### The Test Story

**Client1** (@syftbox.net) has created several folders in their datasite to test different sharing scenarios:
- `private/` - Personal space with user-specific folders
- `public/` - Content everyone can read
- `shared/` - Content shared with specific collaborators
- `dynamic/` - Folder for testing permission changes

The system tests against three clients:
- **client1@syftbox.net**: Primary data owner
- **client2@syftbox.net**: Trusted collaborator
- **bad@syftbox.net**: Potential bad actor (used for security testing)

---

### Test 1: User-Only Permissions (`private/`)

**Goal:** Each user can only access files in their own email-named subfolder

**Permission file:** `datasites/client1@syftbox.net/private/syft.pub.yaml`
```yaml
rules:
  - pattern: '{{.UserEmail}}/*'
    access:
      admin: []
      read:
        - 'USER'
      write:
        - 'USER'
```

**File structure:**
```
datasites/client1@syftbox.net/private/
├── syft.pub.yaml
└── client1@syftbox.net/
    └── secret.txt
```

**Security guarantee:**
- ✅ Client1 accessing `private/client1@syftbox.net/*` → Allowed (template matches their email)
- ❌ Client2 accessing `private/client1@syftbox.net/*` → Denied (template doesn't match)
- ❌ Bad actor accessing `private/client1@syftbox.net/*` → Denied (template doesn't match)

The `{{.UserEmail}}` template ensures user-specific isolation - each user can only access a folder matching their email address.

---

### Test 2: Public Read, Owner Write (`public/`)

**Goal:** Everyone can read files, but only the owner can modify them

**Permission file:** `datasites/client1@syftbox.net/public/syft.pub.yaml`
```yaml
rules:
  - pattern: '**'
    access:
      admin: []
      read:
        - '*'
      write: []
```

**Security guarantee:**
- ✅ Everyone (client1, client2, bad actor) can read files
- ✅ Client1 (owner) can write new files
- ❌ Client2 and bad actor cannot write or modify files

This pattern is ideal for sharing public datasets or results while maintaining control.

---

### Test 3: Specific User Collaboration (`shared/`)

**Goal:** Share with specific trusted users only

**Permission file:** `datasites/client1@syftbox.net/shared/syft.pub.yaml`
```yaml
rules:
  - pattern: '**'
    access:
      admin: []
      read:
        - 'client2@syftbox.net'
      write: []
```

**Security guarantee:**
- ✅ Client1 (owner) has full access (implicit)
- ✅ Client2 can read files
- ❌ Client2 cannot modify files
- ❌ Bad actor cannot read or write

This enables controlled collaboration with specific team members.

---

### Test 4: Dynamic Permission Updates

**Goal:** Test that permission changes propagate correctly

**Initial state (private):**
```yaml
rules:
  - pattern: '{{.UserEmail}}/*'
    access:
      admin: []
      read:
        - 'USER'
      write:
        - 'USER'
```

**Updated state (public):**
```yaml
rules:
  - pattern: '**/*'
    access:
      admin: []
      read:
        - '*'
      write: []
```

**Security guarantee:**
- Permission changes propagate immediately
- Access control is enforced in real-time
- No caching issues that could lead to unauthorized access

---

## YAML Format Requirements

The SyftBox server requires specific YAML formatting with proper indentation:

```yaml
rules:
  - pattern: '**'    # Note the 2-space indentation before the dash
    access:
      admin: []      # Required field, even if empty
      read:
        - '*'        # List items with proper indentation
      write: []      # Empty list means owner-only write
```

**Critical formatting rules:**
1. Use 2 spaces before the `-` in rules array
2. Always include the `admin: []` field
3. Use `**` for current directory and subdirectories
4. Use `*` for current directory only
5. Empty arrays `[]` mean no explicit permissions (owner still has implicit access)

## Key Security Concepts

### Owner Privileges
- Datasite owners always have implicit full access to all their files
- No need to explicitly list the owner in permission rules
- Owner is determined from the first path segment (e.g., `alice` in `/alice/data/`)

### Access Levels
1. **Read**: View file contents and metadata
2. **Write**: Create, modify, or delete files
3. **Admin**: Full control including modifying ACL files

### Special Tokens
- `*`: Wildcard representing all users (public access)
- `USER`: Dynamic token that resolves to the requesting user's email
- `{{.UserEmail}}`: Template that creates user-specific paths

### Pattern Matching
- `*`: Matches files in current directory only
- `**`: Matches files in current directory and all subdirectories
- `**/*`: Matches files in all subdirectories (not current directory)
- `{{.UserEmail}}/*`: Template pattern for user-specific folders

## BioVault-Specific Security Considerations

### Protecting Genomic Data
BioVault will use SyftBox permissions to:
1. Keep raw genomic data private by default
2. Share analysis results selectively with collaborators
3. Publish aggregate statistics publicly while protecting individual data
4. Enable time-limited access for specific research projects

### Compliance and Audit
- All access attempts are logged
- Permission changes are tracked
- Terminal nodes prevent accidental permission escalation
- Regular permission audits can be automated

### Best Practices for BioVault Users

1. **Start Private**: All data should be private by default
2. **Use Terminal Nodes**: For sensitive genomic data directories
3. **Explicit Sharing**: Only share with specific email addresses, avoid wildcards
4. **Regular Audits**: Review permission files periodically
5. **Test Permissions**: Always verify access controls before sharing sensitive data
6. **Template Isolation**: Use `{{.UserEmail}}` patterns for user-specific workspaces

## Testing Security

BioVault includes automated security tests that verify:
- User isolation works correctly
- Public sharing is read-only
- Specific user permissions are enforced
- Bad actors cannot access restricted data
- Permission changes propagate correctly

These tests run nightly in CI to ensure the security model remains intact.

## Reporting Security Issues

If you discover a security vulnerability in BioVault or its use of SyftBox permissions, please report it to:
- Email: security@biovault.org (placeholder)
- Do not disclose security issues publicly until they have been addressed

## Additional Resources

- [SyftBox ACL Documentation](syftbox/docs/acl-system.md)
- [Development Guide](DEV.md)
- [Integration Tests](cli/tests/permission_test.rs)