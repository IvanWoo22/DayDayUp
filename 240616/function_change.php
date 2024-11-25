<?php

function getTestResults($id, $conn): false|string
{
    ob_start();
    ?>
    <div class="container px-3 my-3">
        <?php
        $sql = "SELECT * FROM test_numeric_record WHERE PersonID = ? ORDER BY RecordTime";
        $stmt = $conn->prepare($sql);
        $stmt->bind_param('i', $id);
        $stmt->execute();
        $result = $stmt->get_result();
        if ($result && $result->num_rows > 0) {
            $resultsData = [];
            while ($row = $result->fetch_assoc()) {
                $resultsData[] = $row;
            }

            $groupedData = [];
            foreach ($resultsData as $row) {
                $bigItem = htmlspecialchars($row['BigItem'] ?? '', ENT_QUOTES);
                $item = htmlspecialchars($row['Item'] ?? '', ENT_QUOTES);
                $groupedData[$bigItem][$item][] = $row;
            }


            foreach ($groupedData as $bigItem => $items) {
                if (empty($bigItem)) {
                    continue;
                }
                ?>
                <h3 style="margin:30px"><?= $bigItem ?></h3>
                <div class="accordion" id="<?= $bigItem ?>">
                    <?php
                    foreach ($items as $item => $rows) {
                        if (empty($item)) {
                            continue;
                        }
                        ?>
                        <div class="accordion-item">
                            <h3 class="accordion-header" id="heading<?= $bigItem . "-" . $item ?>">
                                <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse"
                                        data-bs-target="#collapse<?= $bigItem . "-" . $item ?>"
                                        aria-expanded="false"
                                        aria-controls="collapse<?= $bigItem . "-" . $item ?>">
                                    <?= $item ?>
                                </button>
                            </h3>
                            <div id="collapse<?= $bigItem . "-" . $item ?>"
                                 class="accordion-collapse collapse"
                                 aria-labelledby="heading<?= $bigItem . "-" . $item ?>">
                                <div class="accordion-body">
                                    <table class="table table-hover">
                                        <thead>
                                        <tr>
                                            <th scope="col">RecordTime</th>
                                            <th scope="col">Result</th>
                                            <th scope="col">Unit</th>
                                            <th scope="col">RefRange</th>
                                        </tr>
                                        </thead>
                                        <tbody>
                                        <?php
                                        foreach ($rows as $row) {
                                            $recordTime = htmlspecialchars($row['RecordTime'] ?? '', ENT_QUOTES);
                                            $result = htmlspecialchars($row['Result'] ?? '', ENT_QUOTES);
                                            $unit = htmlspecialchars($row['Unit'] ?? '', ENT_QUOTES);
                                            $refRange = htmlspecialchars($row['RefRange'] ?? '', ENT_QUOTES);
                                            ?>
                                            <tr>
                                                <td><?= $recordTime ?></td>
                                                <td><?= $result ?></td>
                                                <td><?= $unit ?></td>
                                                <td><?= $refRange ?></td>
                                            </tr>
                                        <?php } ?>
                                        </tbody>
                                    </table>
                                </div>
                            </div>
                        </div>
                    <?php } ?>
                </div>
                <?php
            }
        } else {
            echo "<p>没有找到相关的化验记录。</p>";
        } ?>
    </div>
    <?php
    return ob_get_clean();
}
