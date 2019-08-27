use erasure_coding::{Encoder, Decoder, Block};

#[test]
fn hello_world() {
    let original_data = b"Hello, World!";
    let encoder = Encoder::new(original_data.len() as u8, 8);
    let decoder = Decoder::new(original_data.len() as u8, 8);

    let (data, repair) = encoder.encode(original_data);

    let mut erased_data: Vec<Option<Block>> = data.into_iter().map(Some).collect();
    erased_data[0] = None;
    erased_data[1] = None;
    erased_data[2] = None;
    erased_data[3] = None;

    let repaired_data = decoder.decode(&erased_data, &repair);
    assert_eq!(repaired_data, original_data);
}

#[test]
fn round_trip() {
    let encoder = Encoder::new(3, 2);
    let decoder = Decoder::new(3, 2);
    let original_data = vec![1, 2, 3];

    let (data, repair) = encoder.encode(&original_data);

    let mut erased_data: Vec<Option<Block>> = data.into_iter().map(Some).collect();
    erased_data[0] = None;
    erased_data[2] = None;

    let repaired_data = decoder.decode(&erased_data, &repair);
    assert_eq!(repaired_data, original_data);
}
