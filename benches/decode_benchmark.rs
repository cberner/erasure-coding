use rand::Rng;
use std::time::Instant;
use erasure_coding::{Encoder, Decoder, Block};

const TARGET_TOTAL_BYTES: usize = 128 * 1024 * 1024;
const DATA_SHARD_COUNTS: [usize; 8] = [3, 7, 10, 30, 70, 100, 150, 200];
const REPAIR_SHARD_COUNTS: [usize; 8] = [2, 2, 3, 5, 10, 10, 15, 20];

fn black_box(value: u64) {
    if value == rand::thread_rng().gen() {
        println!("{}", value);
    }
}

fn benchmark(shard_size: u16) -> u64 {
    let mut black_box_value = 0;
    for (data_shards, repair_shards) in DATA_SHARD_COUNTS.iter().zip(REPAIR_SHARD_COUNTS.iter()) {
        let elements = data_shards * shard_size as usize;
        let mut data: Vec<u8> = vec![0; elements];
        for i in 0..elements {
            data[i] = rand::thread_rng().gen();
        }

        let iterations = TARGET_TOTAL_BYTES / elements;
        let encoder = Encoder::new(*data_shards as u8, *repair_shards as u8);
        let (data_blocks, repair) = encoder.encode(&data);
        let mut erased_datas = vec![];
        for _ in 0..iterations {
            let mut erased: Vec<Option<Block>> = data_blocks.iter().map(|x| Some(x.clone())).collect();
            for _ in 0..*repair_shards {
                let i = rand::thread_rng().gen_range(0, *data_shards);
                erased[i] = None;
            }
            erased_datas.push(erased);
        }
        let now = Instant::now();
        for i in 0..iterations {
            let decoder = Decoder::new(*data_shards as u8, *repair_shards as u8);
            let result = decoder.decode(&erased_datas[i], &repair);
            black_box_value += result[0] as u64;
        }
        let elapsed = now.elapsed();
        let elapsed = elapsed.as_secs() as f64 + elapsed.subsec_millis() as f64 * 0.001;
        let throughput = (elements * iterations) as f64 / 1024.0 / 1024.0 / elapsed;
        println!("data shards = {}, repair shards = {}, decoded {} MB in {:.3}secs, throughput: {:.1}MB/s",
                 data_shards,
                 repair_shards,
                 elements * iterations / 1024 / 1024,
                 elapsed,
                 throughput);
    }

    return black_box_value;
}

fn main() {
    let symbol_size = 512;
    println!("Shard size: {} bytes", symbol_size);
    black_box(benchmark(symbol_size));
}
