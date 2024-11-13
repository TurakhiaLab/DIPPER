import random

def generate_random_sequence(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def mutate_sequence(sequence, mutation_rate):
    mutated = list(sequence)
    for i in range(len(mutated)):
        if random.random() < mutation_rate:
            mutated[i] = random.choice('ACGT')
    return ''.join(mutated)

# Parameters
num_sequences = 50
sequence_length = 5000
mutation_rate = 0.1  # Adjust this to get desired distance

# Generate sequences
ancestral = generate_random_sequence(sequence_length)
sequences = [mutate_sequence(ancestral, mutation_rate) for _ in range(num_sequences)]

# Print sequences
for i, seq in enumerate(sequences):
    print(f">Seq{i+1}")
    print(seq)