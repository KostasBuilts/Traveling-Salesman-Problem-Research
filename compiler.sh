read -r -p "Enter assignment: " user_choice

case "$user_choice" in
    [1])
        source_file="Simulated_Annealing.c"
        ;;
    [2])
        source_file="Quantum_Annealing_Simulation.c"
        ;;
    [3])
        source_file="Hopfield_Network.c"
        ;;
    [4])
        source_file="Genetic_Algorithm.c"
        ;;
    [5])
        source_file="Brute_Force.c"
        ;;
    [6])
        source_file="Ant_Colony_Optimization.c"
        ;;

    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

# Compile the chosen source file
gcc "$source_file" -o program_out -lm

# Check for successful compilation
if [ $? -ne 0 ]; then
    echo "Compilation failed."
    exit 1
fi

# Make executable and run
chmod +x program_out
clear
echo "Running $source_file..."
./program_out