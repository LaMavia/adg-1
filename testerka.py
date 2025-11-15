import sys
import time
import subprocess

def parse_output(file_path) -> dict:
    mapping = dict()
    with open(file_path) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                read_id = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                mapping[read_id] = (start, end)
    return mapping

def evaluate_mapping(predicted, truth, tolerance=20) -> dict:
    total = len(truth)
    correct = 0
    incorrect = 0
    unmapped = 0

    for read, (true_start, true_end) in truth.items():
        if read not in predicted:
            unmapped += 1
        else:
            pred_start, pred_end = predicted[read]
            if abs(pred_start - true_start) <= tolerance and abs(pred_end - true_end) <= tolerance:
                correct += 1
            else:
                incorrect += 1

    mapped = total - unmapped

    return {
        "total_reads": total,
        "mapped": mapped,
        "unmapped": unmapped,
        "correct": correct,
        "incorrect": incorrect,
        "mapping_percent": mapped / total * 100,
        "correct_percent": correct / total * 100,
        "incorrect_percent": incorrect / total * 100,
    }


def main(): 
    if len(sys.argv) != 5:
        print("usage: python3 testerka.py reference.fasta reads.fasta ground_truth.txt mapper.py")
        sys.exit(1)

    reference_file, reads, ground_truth, mapper = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    output_file = "output_.txt"

    # run the mapper 
    start_time = time.time()
    subprocess.run(["python3", mapper, reference_file, reads, output_file], check=True, stderr=sys.stderr)
    end_time = time.time() - start_time 

    # get the real and predicted reads
    predicted = parse_output(output_file)
    truth = parse_output(ground_truth)

    results = evaluate_mapping(predicted, truth)
    total_reads = results["total_reads"]
    avg_time = end_time / total_reads if total_reads > 0 else 0

    for key, result in results.items(): 
        res = f"{key}: {result}"
        if "percent" in key:
            res += "%"
        print(res)

    print(f"\nTotal time: {end_time:.2f} s")
    print(f"Average time per read: {avg_time:.3f} s/read")


if __name__ == "__main__":
    main()
