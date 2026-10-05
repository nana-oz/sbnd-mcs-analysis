#!/bin/bash

# 1. ターミナルの Python エイリアスを解除（Conda の Python を優先させる）
unalias python 2>/dev/null

# 2. 古い PYTHONPATH の設定をリセット
unset PYTHONPATH

# 3. Dedicated の Conda 環境を有効化
conda activate sbnd_env

# 4. プロジェクトルート（現在地）を Python の参照パスに追加
export PYTHONPATH="$PWD:$PYTHONPATH"

# 5. 起動確認ログの出力
echo "=========================================="
echo "  sbnd_env setup complete"
echo "  Python: $(which python)"
echo "  PYTHONPATH: $PYTHONPATH"
echo "=========================================="