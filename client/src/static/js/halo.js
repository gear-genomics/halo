import axios from 'axios'
import pako from 'pako'

const urlExample = document.getElementById('link-example').href
var exampleData

$('#mainTab a').on('click', function(e) {
  e.preventDefault()
  $(this).tab('show')
})

const fileInput = document.getElementById('inputFile')
const resultLink = document.getElementById('link-results')
const resultInfo = document.getElementById('result-info')
const sampleSelect = document.getElementById('selectSample')
const resultContainer = document.getElementById('result-container')
const chartsContainer = document.getElementById('charts-container')

const submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', function() {
  resultLink.click()
  run()
})

const exampleButton = document.getElementById('btn-example')
exampleButton.addEventListener('click', showExample)

sampleSelect.addEventListener('change', handleChangeSample)

var inputData

function run() {
  const file = fileInput.files[0]

  // TODO error handling
  if (file === undefined) return

  hideElement(resultContainer)
  chartsContainer.innerHTML = ''
  showElement(resultInfo)

  const fileReader = new FileReader()
  const isGzip = /\.gz$/.test(file.name)
  if (isGzip) {
    fileReader.readAsArrayBuffer(file)
  } else {
    fileReader.readAsText(file)
  }
  fileReader.onload = event => {
    let content = event.target.result
    if (isGzip) {
      content = pako.ungzip(content, { to: 'string' })
    }
    inputData = JSON.parse(content).data
    setupSampleSelect()
    setTimeout(() => {
      vis(inputData[0])
    }, 0)
  }
}

function setupSampleSelect() {
  const samples = inputData.map(rec => rec.sample)
  sampleSelect.innerHTML = samples
    .map((sample, index) => `<option value="${index}">${sample}</option>`)
    .join('\n')
}

function handleChangeSample(event) {
  const index = Number.parseInt(event.target.value)
  hideElement(resultContainer)
  chartsContainer.innerHTML = ''
  showElement(resultInfo)
  setTimeout(() => {
    vis(inputData[index])
  }, 0)
}

function vis(data) {
  showElement(resultContainer)
  for (const coverage of data.coverages) {
    const coverageContainer = document.createElement('div')
    const chartContainer = document.createElement('div')
    coverageContainer.appendChild(chartContainer)
    chartsContainer.appendChild(coverageContainer)
    const title = `${coverage.chromosome} – ${data.sample}`
    renderStrandedCoverageChart(chartContainer, coverage, title)
  }
  hideElement(resultInfo)
}

function renderStrandedCoverageChart(container, data, title) {
  const trace1 = {
    x: data.positions,
    y: data.counts[0].values,
    name: data.counts[0].label,
    type: 'bar'
  }

  const trace2 = {
    x: data.positions,
    y: data.counts[1].values,
    yaxis: 'y2',
    name: data.counts[1].label,
    type: 'bar'
  }

  const traces = [trace1, trace2]

  const layout = {
    title: title || '',
    yaxis: {
      title: 'count',
      domain: [0.54, 1]
    },
    yaxis2: {
      domain: [0, 0.46],
      autorange: 'reversed'
    }
  }

  Plotly.newPlot(container, traces, layout)
}

function showExample() {
  hideElement(resultContainer)
  chartsContainer.innerHTML = ''
  showElement(resultInfo)
  resultLink.click()

  if (exampleData) {
    inputData = exampleData
    setupSampleSelect()
    setTimeout(() => {
      vis(inputData[0])
    }, 0)
  } else {
    axios
      .get(urlExample, {
        responseType: 'arraybuffer'
      })
      .then(response => {
        const content = pako.ungzip(response.data, { to: 'string' })
        exampleData = JSON.parse(content).data
        inputData = exampleData
        setupSampleSelect()
        setTimeout(() => {
          vis(inputData[0])
        }, 0)
      })
      .catch(error => {
        // FIXME proper error handling
        console.error(error)
      })
  }

}

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}
