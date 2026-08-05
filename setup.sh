#!/bin/bash
git init

random_hex() {
    printf "%02x" $((RANDOM % 256))
}

echo "Generating GitHub timeline for 2025-2026..."
echo "This will take about 30-60 seconds..."

for i in {1..2500}; do
    year=$((2025 + RANDOM % 2))
    
    if [ $year -eq 2026 ]; then
        month=$((1 + RANDOM % 8))
        if [ $month -eq 8 ]; then
            day=$((1 + RANDOM % 5))
        else
            day=$((1 + RANDOM % 28))
        fi
    else
        month=$((1 + RANDOM % 12))
        day=$((1 + RANDOM % 28))
    fi
    
    hour=$(printf "%02d" $((RANDOM % 24)))
    minute=$(printf "%02d" $((RANDOM % 60)))
    second=$(printf "%02d" $((RANDOM % 60)))
    
    date="$year-$month-${day}T$hour:$minute:$second.000Z"
    
    commits=$((1 + RANDOM % 4))
    
    for ((j=0; j<commits; j++)); do
        hex=$(random_hex)
        echo -n "$hex" >> dump.txt
        git add dump.txt
        GIT_COMMITTER_DATE="$date" git commit -m "add: $hex" --date="$date" --allow-empty-message > /dev/null 2>&1
    done
done

echo "✅ Timeline generated!"
echo "Total commits: $(git rev-list --count HEAD)"
