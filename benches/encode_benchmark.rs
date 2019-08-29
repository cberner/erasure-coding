use rand::Rng;
use std::time::Instant;
use erasure_coding::Encoder;

const TARGET_TOTAL_BYTES: usize = 128 * 1024 * 1024;
const DATA_SHARD_COUNTS: [usize; 8] = [3, 7, 10, 30, 70, 100, 150, 200];
const REPAIR_SHARD_COUNTS: [usize; 8] = [2, 2, 3, 5, 10, 10, 15, 20];

fn black_box(value: u64) {
    if value == rand::thread_rng().gen() {
        println!("{}", value);
    }
}

fn main() {
    let mut black_box_value = 0;
    let shard_size = 512;
    println!("Shard size: {} bytes", shard_size);
    for (data_shards, repair_shards) in DATA_SHARD_COUNTS.iter().zip(REPAIR_SHARD_COUNTS.iter()) {
        let elements = data_shards * shard_size as usize;
        let mut data: Vec<u8> = vec![0; elements];
        for i in 0..elements {
            data[i] = rand::thread_rng().gen();
        }

        let now = Instant::now();
        let iterations = TARGET_TOTAL_BYTES / elements;
        let encoder = Encoder::new(*data_shards as u8, *repair_shards as u8);
        for _ in 0..iterations {
            let (_, repair) = encoder.encode(&data);
            black_box_value += repair[0][0] as u64;
        }
        let elapsed = now.elapsed();
        let elapsed = elapsed.as_secs() as f64 + elapsed.subsec_millis() as f64 * 0.001;
        let throughput = (elements * iterations) as f64 / 1024.0 / 1024.0 / elapsed;
        println!("data shards = {}, repair shards = {}, encoded {} MB in {:.3}secs, throughput: {:.1}MB/s",
                 data_shards,
                 repair_shards,
                 elements * iterations / 1024 / 1024,
                 elapsed,
                 throughput);
    }
    black_box(black_box_value);
}
