VSCODE_DIR=/var/tmp/pelai/vscode-server/.vscode-server
COMMIT=8b640eef5a6c6089c029249d48efa5c99adf7d51

# 先殺掉殘留 VS Code server
pkill -u $USER -f "$VSCODE_DIR" || true
pkill -u $USER -f code-server || true
pkill -u $USER -f vscode-server || true

# 清掉壞掉的 server 與 staging
rm -rf "$VSCODE_DIR/cli/servers/Stable-$COMMIT"
rm -rf "$VSCODE_DIR/cli/servers/Stable-$COMMIT.staging"

# 也清掉可能殘留的 tmp
rm -rf /tmp/.tmp*/vscode-server-linux-x64.tar.gz 2>/dev/null
rm -rf /tmp/code-* 2>/dev/null
