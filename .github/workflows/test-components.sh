#!/bin/bash
set -e

echo "=== Testing Individual Components ==="

# Test 1: Check if mfc.sh exists and is executable
echo "1. Testing mfc.sh..."
if [ -f "./mfc.sh" ]; then
    chmod +x ./mfc.sh
    echo "✓ mfc.sh found and made executable"
else
    echo "✗ mfc.sh not found - creating mock for testing"
    cat > ./mfc.sh << 'MOCK'
#!/bin/bash
echo "Mock benchmark running..."
mkdir -p bench
echo "name: test_benchmark" > bench/test.yaml
echo "time: 1.23" >> bench/test.yaml
echo "Mock benchmark completed"
MOCK
    chmod +x ./mfc.sh
fi

# Test 2: Run benchmark command
echo "2. Testing benchmark execution..."
./mfc.sh bench -o bench || echo "Benchmark failed, continuing with mock data"

# Test 3: Check yq installation and YAML to JSON conversion
echo "3. Testing YAML to JSON conversion..."
if ! command -v yq &> /dev/null; then
    echo "Installing yq..."
    sudo wget -qO /usr/local/bin/yq https://github.com/mikefarah/yq/releases/latest/download/yq_linux_amd64
    sudo chmod +x /usr/local/bin/yq
fi

# Convert any YAML files found
find . -maxdepth 1 -name "*.yaml" -exec sh -c 'echo "Converting $1"; yq eval -o=json "$1" > "${1%.yaml}.json"' _ {} \;

# Test 4: Check JSON output format
echo "4. Validating JSON format..."
for json_file in *.json; do
    if [ -f "$json_file" ]; then
        echo "Checking $json_file:"
        python3 -c "import json; json.load(open('$json_file')); print('✓ Valid JSON')" || echo "✗ Invalid JSON"
    fi
done

echo "=== Component Testing Complete ==="
