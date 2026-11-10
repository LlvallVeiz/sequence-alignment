.PHONY: build test vet fmt examples clean

build:
	go build -o seqalign ./cmd/seqalign

test:
	go test ./...

vet:
	go vet ./...

fmt:
	gofmt -l .

# Regenerate the example outputs in examples/output/ from real runs.
examples: build
	sh examples/run_examples.sh

clean:
	rm -f seqalign
