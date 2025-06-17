command=$@
  
$command &
pid=$!

# Initialize VmPeak variable
vmpeak=0

# Poll the process status until it finishes
while [ -e /proc/$pid ]; do
    current_vmpeak=$(grep VmPeak /proc/$pid/status | awk '{print $2}')
    if [ -n "$current_vmpeak" ] && [ "$current_vmpeak" -gt "$vmpeak" ]; then
        vmpeak=$current_vmpeak
    fi
    sleep 0.1
done

printf "VmPeak:\t%s kB\n" "$vmpeak" >&2