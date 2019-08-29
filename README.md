[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Fcberner%2Ferasure-coding.svg?type=shield)](https://app.fossa.io/projects/git%2Bgithub.com%2Fcberner%2Ferasure-coding?ref=badge_shield)

## Benchmarks

The following were run on an Intel Core i5-6600K @ 3.50GHz

```
Shard size: 512 bytes
data shards = 3, repair shards = 2, decoded 127 MB in 0.206secs, throughput: 621.4MB/s
data shards = 7, repair shards = 2, decoded 127 MB in 0.141secs, throughput: 907.8MB/s
data shards = 10, repair shards = 3, decoded 127 MB in 0.149secs, throughput: 859.0MB/s
data shards = 30, repair shards = 5, decoded 127 MB in 0.138secs, throughput: 927.5MB/s
data shards = 70, repair shards = 10, decoded 127 MB in 0.187secs, throughput: 684.3MB/s
data shards = 100, repair shards = 10, decoded 127 MB in 0.173secs, throughput: 739.8MB/s
data shards = 150, repair shards = 15, decoded 127 MB in 0.230secs, throughput: 556.3MB/s
data shards = 200, repair shards = 20, decoded 127 MB in 0.271secs, throughput: 472.1MB/s

Shard size: 512 bytes
data shards = 3, repair shards = 2, encoded 127 MB in 0.056secs, throughput: 2285.7MB/s
data shards = 7, repair shards = 2, encoded 127 MB in 0.044secs, throughput: 2909.1MB/s
data shards = 10, repair shards = 3, encoded 127 MB in 0.047secs, throughput: 2723.4MB/s
data shards = 30, repair shards = 5, encoded 127 MB in 0.052secs, throughput: 2461.5MB/s
data shards = 70, repair shards = 10, encoded 127 MB in 0.077secs, throughput: 1661.9MB/s
data shards = 100, repair shards = 10, encoded 127 MB in 0.076secs, throughput: 1683.9MB/s
data shards = 150, repair shards = 15, encoded 127 MB in 0.098secs, throughput: 1305.7MB/s
data shards = 200, repair shards = 20, encoded 127 MB in 0.131secs, throughput: 976.6MB/s
```
## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.


[![FOSSA Status](https://app.fossa.io/api/projects/git%2Bgithub.com%2Fcberner%2Ferasure-coding.svg?type=large)](https://app.fossa.io/projects/git%2Bgithub.com%2Fcberner%2Ferasure-coding?ref=badge_large)

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.