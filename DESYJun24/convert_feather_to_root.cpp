#include <iostream>
#include <vector>
#include <map>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <arrow/ipc/api.h>

arrow::Status ConvertFeatherToRoot(const std::string& featherFile, const std::string& rootFile) {
    // Open Feather file
    std::shared_ptr<arrow::io::ReadableFile> infile;
    ARROW_ASSIGN_OR_RAISE(infile, arrow::io::ReadableFile::Open(featherFile));

    std::shared_ptr<arrow::ipc::RecordBatchFileReader> reader;
    ARROW_ASSIGN_OR_RAISE(reader, arrow::ipc::RecordBatchFileReader::Open(infile));

    std::vector<std::shared_ptr<arrow::RecordBatch>> batches;
    for (int i = 0; i < reader->num_record_batches(); i++) {
        ARROW_ASSIGN_OR_RAISE(auto batch, reader->ReadRecordBatch(i));
        batches.push_back(batch);
    }

    // Combine batches into a table
    std::shared_ptr<arrow::Table> table;
    ARROW_ASSIGN_OR_RAISE(table, arrow::Table::FromRecordBatches(batches));

    // Extract column indices
    std::shared_ptr<arrow::ChunkedArray> evt_column = table->GetColumnByName("evt");
    if (!evt_column) {
        std::cerr << "Column 'evt' not found!" << std::endl;
        return arrow::Status::OK();
    }

    // Open ROOT file and create TTree
    TFile f(rootFile.c_str(), "RECREATE");
    TTree tree("Events", "Tree storing feather data");

    std::map<int, std::vector<int>> evt_data;
    std::map<std::string, std::map<int, std::vector<int>>> int_columns;
    std::map<std::string, std::map<int, std::vector<float>>> float_columns;

    // Column names
    std::vector<std::string> int_col_names = {"bcid", "l1a_counter", "ea", "board", "row", "col", "toa", "tot", "cal"};
    std::vector<std::string> float_col_names = {"acal"};

    // Read columns into maps
    for (const auto& name : int_col_names) {
        int_columns[name] = {};
        auto col = table->GetColumnByName(name);
        if (!col) continue;
        auto int_array = std::static_pointer_cast<arrow::Int64Array>(col->chunk(0));
        for (int i = 0; i < int_array->length(); i++) {
            int evt = evt_column->chunk(0)->GetScalar(i).ValueOrDie()->hash();
            int val = int_array->Value(i);
            int_columns[name][evt].push_back(val);
        }
    }

    for (const auto& name : float_col_names) {
        float_columns[name] = {};
        auto col = table->GetColumnByName(name);
        if (!col) continue;
        auto float_array = std::static_pointer_cast<arrow::DoubleArray>(col->chunk(0));
        for (int i = 0; i < float_array->length(); i++) {
            int evt = evt_column->chunk(0)->GetScalar(i).ValueOrDie()->hash();
            float val = static_cast<float>(float_array->Value(i));
            float_columns[name][evt].push_back(val);
        }
    }

    // Create branches in the tree
    std::vector<int> evt_vector;
    std::map<std::string, std::vector<int>> int_vectors;
    std::map<std::string, std::vector<float>> float_vectors;

    tree.Branch("evt", &evt_vector);
    for (const auto& name : int_col_names) {
        tree.Branch(name.c_str(), &int_vectors[name]);
    }
    for (const auto& name : float_col_names) {
        tree.Branch(name.c_str(), &float_vectors[name]);
    }

    // Fill the tree
    for (const auto& [evt, _] : evt_data) {
        evt_vector = {evt};
        for (const auto& name : int_col_names) {
            int_vectors[name] = int_columns[name][evt];
        }
        for (const auto& name : float_col_names) {
            float_vectors[name] = float_columns[name][evt];
        }
        tree.Fill();
    }

    tree.Write();
    f.Close();
    std::cout << "ROOT file created successfully: " << rootFile << std::endl;

    return arrow::Status::OK();
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
	std::cout << "You must enter two arguments" << std::endl; 
	return 1;
    }
    std::string inFile = argv[1];
    std::string rootFile = argv[2];
    arrow::Status status = ConvertFeatherToRoot(inFile, rootFile);
    if (!status.ok()) {
        std::cerr << "Error: " << status.message() << std::endl;
        return 1;
    }
    return 0;
}
